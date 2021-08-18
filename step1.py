import argparse
import requests
import logging

from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA

import pandas as pd
import numpy as np
from difflib import SequenceMatcher

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def fetch_ensembl_transcript(ensembl_transcript_id, exon_annot = False):
    """Fetch the requested Ensembl transcript.
    Get the requested Ensembl transcript, together with exon and
    coding region (CDS) boundaries.
    Parameters
    ----------
    ensembl_transcript_id : str
      the ensembl transcript id, of the form ENST...

    Returns
    -------
    `Bio.SeqRecord`
      The requested transcript sequence, in 5' -> 3' order, together
      with exon and CDS features. The coordinates of exons and CDS
      features are relative to the sequence fragment.
    """

    # TODO: Validate ensembl_transcript_id is a valid transcript id

    base_url = "http://rest.ensembl.org"

    # First, fetch the transcript sequence
    url = base_url + f"/sequence/id/{ensembl_transcript_id}"

    log.info(f"Querying Ensembl for sequence of {ensembl_transcript_id}")
    response = requests.get(url, { "type": "genomic",
                                   "content-type": "application/json" })

    try:
        response.raise_for_status()
    except requests.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
        raise ValueError(reponse.text)

    response_data = response.json()

    try:
        description = response_data['desc'].split(':')
        species = description[1]
        try:
            chromosome_number = int(description[2])
        except:
            chromosome_number = str(description[2])
        sequence_left = int(description[3])
        sequence_right = int(description[4])
        transcript_strand = int(description[5])

        if sequence_left > sequence_right:
            raise ValueError(f"Expected left sequence boundary {sequence_left} "
                             f"<= right sequence boundary {sequence_right}: did "
                             "the format of the Ensembl REST response change?")

        sequence_id = response_data['id']

        seq_str = response_data['seq']

        log.info(f"Retrieved sequence {response_data['desc']} of length "
                 f"{sequence_right - sequence_left} for species {species} on "
                 f"strand {transcript_strand}")
    except (KeyError, ValueError) as e:
        log.error(e)
        log.error('Error parsing sequence metadata from Ensembl REST response - '
                  'did the format of the response change?')
        raise ValueError(e)

    seq = Seq(seq_str, IUPACUnambiguousDNA())

    record = SeqRecord(seq, id=sequence_id,
                       description=":".join(description))
    if exon_annot:

        url = base_url + f"/overlap/id/{ensembl_transcript_id}"

        log.info(f"Querying Ensembl for overlaps of {ensembl_transcript_id}")
        response = requests.get(url, { "feature": ["cds", "exon"],
                                       "content-type": "application/json" })
        try:
            response.raise_for_status()
        except requests.HTTPError:
            log.error("Ensembl sequence REST query returned error "
                      "{}".format(response.text))
            raise ValueError(reponse.text)

        response_data = response.json()

        try:
            # Handle the unlikely event of a single piece of information
            # overlapping a lonely transcript
            if not hasattr(response_data, '__iter__'):
                response_data = [response_data]

            for response_datum in response_data:
                if response_datum['Parent'] != ensembl_transcript_id:
                    continue

                if response_datum['assembly_name'] != species:
                    continue

                # We store feature locations 0-indexed from the left-most
                # sequence boundary
                record.features.append(SeqFeature(
                    location=FeatureLocation(
                        int(response_datum['start']) - sequence_left,
                        int(response_datum['end']) - sequence_left + 1,
                        strand=int(response_datum['strand'])),
                    type=response_datum['feature_type']))
            num_exon_boundaries = len([f for f in record.features
                                       if f.type == 'exon'])

            num_cds_boundaries = len([f for f in record.features
                                      if f.type == 'cds'])

            log.info(f"Retrieved {num_exon_boundaries} exons and "
                     f"{num_cds_boundaries} coding regions for transcript "
                     f"{ensembl_transcript_id}")
        except (KeyError, ValueError) as e:
            log.error(e)
            log.error('Error parsing overlap metadata from Ensembl REST response - '
                      'did the format of the response change?')
            raise ValueError(e)

    record.annotations['reference_species'] = species
    record.annotations['reference_chromosome_number'] = chromosome_number
    record.annotations['reference_left_index'] = sequence_left
    record.annotations['reference_right_index'] = sequence_right
    record.annotations['transcript_strand'] = transcript_strand

    # Finally, sort features by their start locations
    record.features.sort(key=lambda f: f.location.start)

    return record

def fetch_gggenome_match(seq, genome='hg38', mismatches=0):
    """Returns matches from gggenome service. For example:
    Parameters
    ----------
    genome : str
        See https://gggenome.dbcls.jp/help.html
    mismatches : int
        Number of allowed mismatches or gaps. IMPORTANT:
        `mismatches=4` causes timeouts on 20bp and
        `mismatches=2` takes ~30s per 20bp request in testing.
    strand : str
        '+', '-', or '', indicating whether to search the reference
        strand, the anti-reference strand, or both strands.
    Returns
    -------
    `List[Dict]`
        For example: ```
            [{
                "align": "|||||| |||||| || || |||",
                "del": 0,
                "edit": "======I======I==X==X===",
                "ins": 2,
                "match": 19,
                "mis": 2,
                "name": "chr1",
                "position": 161696498,
                "position_end": 161696518,
                "query": "TTCCGGCGCGCCGAGTCCTTAGG",
                "sbjct": "TTCCGG-GCGCCG-GTGCTGAGG",
                "snippet": "TGTGGGTCTGCAAGGAGCCCTCGCGGGAAGCAGGAAGGAGCGGGGTCGCGGAGCGGTGGACAAGCCGGCGCCGTTGCTCCCCGCCCTCTCCGTAGAGCTGTTCCGGGCGCCGGTGCTGAGGGTGATGGGTCCGCGGGAGGCCCGCGGCGCGGCGCTGGGTGGGGTGGTGCTGCGCTGCGACACGCGCCTGCACCCGCAGAAGCGCGACACGCCGCTGCAGT",
                "snippet_end": 161696618,
                "snippet_pos": 161696398,
                "strand": "+"
            }]
    ```
    """

    url = 'http://gggenome.dbcls.jp/{}/{}/{}.json'.format(
        genome, mismatches, seq)

    log.info('Querying gggenome for offtargets')
    response = requests.get(url)

    try:
        response.raise_for_status()
    except requests.HTTPError:
        log.error(f"gggenome REST query returned error {response.text}")
        raise ValueError(response.text)

    data = response.json()

    if data['error'] != 'none':
        raise RuntimeError('gggenome error: "{}"'.format(data['error']))

    return data['results']

def fetch_ensembl_sequence(chromosome, region_left, region_right, expand = 200):
    '''
    Returns genome sequence based on chromosome range. The sequence is expanded by a flat amount on both the 5 and
    3 termini.
    '''
    base_url = "http://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{region_left}..{region_right}:1?expand_5prime={expand};expand_3prime={expand}"
    r = requests.get(base_url + ext, headers = {"Content-Type":"text/plain"})

    if not r.ok:
        r.raise_for_status()

    sequence = Seq(r.text, IUPACUnambiguousDNA())
    return sequence

def soft_match(ultramer, expanded_sequence):
    '''
    In situation where ultramer that aligns with expanded_sequence perfectly on one side, and
    is expected to align with small mismatches on the other side, this function returns
    the ultramer boundaries
    '''
    if expanded_sequence.find(ultramer[:16]) + 1:
        leftend = expanded_sequence.find(ultramer[:16])
        left_matches = SequenceMatcher(None, ultramer[-16:], expanded_sequence[leftend + 55:leftend+120]).get_matching_blocks()
    elif expanded_sequence.find(ultramer.reverse_complement()[:16]) + 1:
        leftend = expanded_sequence.find(ultramer.reverse_complement()[:16])
        left_matches = SequenceMatcher(None, ultramer.reverse_complement()[-16:], expanded_sequence[leftend + 55:leftend+120]).get_matching_blocks()

    try:
        match_highest = left_matches[0]
        for match in left_matches:
            if match.size > 0 and match.b > match_highest.b:
                match_highest = match
        assert match_highest.a + match_highest.size == 16, 'soft_match hasn\'t aligned whole ultramer segment'
        return [leftend, leftend + 55 + match_highest.b + match_highest.size]
    except:
        pass

    if expanded_sequence.find(ultramer[-16:]) + 1:
        rightend = expanded_sequence.find(ultramer[-16:]) + 16
        right_matches = SequenceMatcher(None, ultramer[:16], expanded_sequence[rightend - 120:rightend - 55]).get_matching_blocks()
    elif expanded_sequence.find(ultramer.reverse_complement()[-16:]) + 1:
        rightend = expanded_sequence.find(ultramer.reverse_complement()[-16:]) + 16
        right_matches = SequenceMatcher(None, ultramer.reverse_complement()[:16], expanded_sequence[rightend - 120:rightend - 55]).get_matching_blocks()

    try:
        match_lowest = right_matches[0]
        match_highest = right_matches[0]
        for match in right_matches:
            if match.size > 0 and match.b > match_highest.b:
                match_highest = match
        assert match_highest.a + match_highest.size == 16, 'soft_match hasn\'t aligned whole ultramer segment'
        assert match_lowest.size > 0, 'check match_lowest'
        return [rightend - 120 + match_lowest.b, rightend]
    except:
        pass

def delimit_insertion(platepath):
    '''
    This converts protospacers and ultramers designed against gene transcripts into gene coordinate,
    which are to be fed into crispr-primer for primer design.
    '''
    ultramersdf = pd.read_excel(platepath, index = None)

    assert {'transcript','gene','protospacer','Ultramer'}.issubset(set(ultramersdf.columns)),'excel header columns should include "transcript","gene","protospacer",and "Ultramer"'


    platedf = pd.DataFrame(columns=['sample','chromosome','ultramer_range_left','ultramer_range_right','bed_range'])

    for index, row in ultramersdf.iterrows():
        if row['transcript'] == row['transcript']:
            expand = 200

            transcript = row['transcript']
            transcript = transcript.split()[0]
            assert transcript[:4]=='ENST' and len(transcript) == 15, 'check transcript ID formatting'

            well = row['well']

            protospacer = row['protospacer']
            protospacer = Seq(protospacer.upper(),IUPACUnambiguousDNA())

            ultramer = row['Ultramer']
            ultramer = Seq(ultramer.upper(), IUPACUnambiguousDNA())

            record = fetch_ensembl_transcript(transcript)

            chromosome, region_left, region_right = (record.annotations['reference_chromosome_number'],
            record.annotations['reference_left_index'], record.annotations['reference_right_index'])

            sequence = record.seq
            expanded_sequence = fetch_ensembl_sequence(chromosome, region_left, region_right, expand)

            #Note that expanded_sequence will always be in the direction of the reference genome, and has no
            #bearing on the strandedness of the transcript. To retrieve that information, use
            #record.annotations['transcript strand']

            assert (protospacer in expanded_sequence) or (protospacer.reverse_complement() in expanded_sequence), f'{index} protospacer not found in transcript'


            if -1 not in [expanded_sequence.find(ultramer[:25]), expanded_sequence.find(ultramer[-25:])]:
                ult_range = [expanded_sequence.find(ultramer[:25]), expanded_sequence.find(ultramer[-25:])+25]
            elif -1 not in [expanded_sequence.find(ultramer.reverse_complement()[:25]),
                             expanded_sequence.find(ultramer.reverse_complement()[-25:])]:
                ult_range = [expanded_sequence.find(ultramer.reverse_complement()[:25]),
                             expanded_sequence.find(ultramer.reverse_complement()[-25:])+25]
            else:
                print(well, 'can\'t find 25bp ends of ultramer, trying 16')
                if -1 not in [expanded_sequence.find(ultramer[:16]), expanded_sequence.find(ultramer[-16:])]:
                    ult_range = [expanded_sequence.find(ultramer[:16]), expanded_sequence.find(ultramer[-16:])+16]
                elif -1 not in [expanded_sequence.find(ultramer.reverse_complement()[:16]),
                                 expanded_sequence.find(ultramer.reverse_complement()[-16:])]:
                    ult_range = [expanded_sequence.find(ultramer.reverse_complement()[:16]),
                                 expanded_sequence.find(ultramer.reverse_complement()[-16:])+16]
                else:
                    print(well, 'can\'t find 16bp ends of ultramer, looking for soft alignment')
                    ult_range = soft_match(ultramer, expanded_sequence)





            assert ult_range[1] - ult_range[0] in range(70,140), 'did we change the total length of homology arms?'

            check_strand_consistency(well, expanded_sequence, protospacer, ultramer)

            ultramer_range_left = region_left - expand + ult_range[0]
            ultramer_range_right = region_left - expand + ult_range[1]

            rowdf = pd.DataFrame([[well, chromosome, ultramer_range_left, ultramer_range_right, f'chr{chromosome}:{ultramer_range_left}-{ultramer_range_right}']],
                                columns = ['sample','chromosome','ultramer_range_left','ultramer_range_right','bed_range'])
            platedf = platedf.append(rowdf)

        #added to search for Jin protospacers
        elif row['protospacer'] == row['protospacer']:
            well = row['well']

            expand = 500

            protospacer = row['protospacer']

            query_results = fetch_gggenome_match(protospacer)

            protospacer = Seq(protospacer.upper(),IUPACUnambiguousDNA())
            ultramer = row['Ultramer']
            ultramer = Seq(ultramer.upper(), IUPACUnambiguousDNA())

            ult_range = []
            for query_result in query_results:
                if not ult_range:
                    chromosome, region_left, region_right = query_result['name'], query_result['position'], query_result['position_end']
                    assert 'chr' in chromosome, 'check gggenome output'

                    expanded_sequence = fetch_ensembl_sequence(chromosome, region_left, region_right, expand)

                    assert (protospacer in expanded_sequence) or (protospacer.reverse_complement() in expanded_sequence), f'{index} protospacer not found in transcript'

                    if -1 not in [expanded_sequence.find(ultramer[:25]), expanded_sequence.find(ultramer[-25:])]:
                        ult_range = [expanded_sequence.find(ultramer[:25]), expanded_sequence.find(ultramer[-25:])+25]
                    elif -1 not in [expanded_sequence.find(ultramer.reverse_complement()[:25]),
                                     expanded_sequence.find(ultramer.reverse_complement()[-25:])]:
                        ult_range = [expanded_sequence.find(ultramer.reverse_complement()[:25]),
                                     expanded_sequence.find(ultramer.reverse_complement()[-25:])+25]
                    else:
                        print(well, 'can\'t find 25bp ends of ultramer, trying 16')
                        if -1 not in [expanded_sequence.find(ultramer[:16]), expanded_sequence.find(ultramer[-16:])]:
                            ult_range = [expanded_sequence.find(ultramer[:16]), expanded_sequence.find(ultramer[-16:])+16]
                        elif -1 not in [expanded_sequence.find(ultramer.reverse_complement()[:16]),
                                         expanded_sequence.find(ultramer.reverse_complement()[-16:])]:
                            ult_range = [expanded_sequence.find(ultramer.reverse_complement()[:16]),
                                         expanded_sequence.find(ultramer.reverse_complement()[-16:])+16]
                        else:
                            print(well, 'can\'t find ultramer')

            assert ult_range[1] - ult_range[0] in range(70,150), 'did we change the total length of homology arms?'

            check_strand_consistency(well, expanded_sequence, protospacer, ultramer)

            ultramer_range_left = region_left - expand + ult_range[0]
            ultramer_range_right = region_left - expand + ult_range[1]

            chromosome = chromosome[3:]

            rowdf = pd.DataFrame([[well, chromosome, ultramer_range_left, ultramer_range_right, f'chr{chromosome}:{ultramer_range_left}-{ultramer_range_right}']],
                                columns = ['sample','chromosome','ultramer_range_left','ultramer_range_right','bed_range'])
            platedf = platedf.append(rowdf)

    return platedf


def check_strand_consistency(well, expanded_sequence, protospacer, ultramer):
    '''
    will print note if protospacer not on same strand as ultramer
    '''
    if protospacer in expanded_sequence:
        break_strand = 1
    else:
        break_strand = 0

    if [expanded_sequence.find(ultramer[:25]), expanded_sequence.find(ultramer[-25:])] != [-1, -1]:
        ult_strand = 1
    elif [expanded_sequence.find(ultramer.reverse_complement()[:25]),
                     expanded_sequence.find(ultramer.reverse_complement()[-25:])] != [-1, -1]:
        ult_strand = 0
    elif [expanded_sequence.find(ultramer[:16]), expanded_sequence.find(ultramer[-16:])] != [-1, -1]:
        ult_strand = 1
    elif [expanded_sequence.find(ultramer.reverse_complement()[:16]),
                     expanded_sequence.find(ultramer.reverse_complement()[-16:])] != [-1, -1]:
        ult_strand = 0
    if break_strand != ult_strand:
        print(well, f"strandedness inconsistent")

if __name__ == "__main__":
    platedf = delimit_insertion(sys.argv[0])
    platedf.to_csv('mNGplate_primers_in.csv',columns=['sample', 'bed_range'], header = None, index=None)
    
