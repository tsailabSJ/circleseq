from __future__ import print_function

import subprocess
import sys
import os
import argparse
import regex
import re
import HTSeq
import pyfaidx
from findCleavageSites import get_sequence, regexFromSequence, alignSequences, reverseComplement, extendedPattern, realignedSequences

"""
Run samtools:mpileup and get all identified variants in the window sequences
"""
def snpCall(basename, matched_file, reference_genome, bam_file, out_snp, search_radius):
    output_folder = os.path.dirname(out_snp)

    # open matched file
    regions = list()
    with open(matched_file, 'rU') as f:
        f.readline()
        for line in f:
            site = line.strip().split('\t')
            regions.append([site[0], int(site[1]) - search_radius, int(site[2]) + search_radius, '*', bam_file, ':'.join([site[24], site[3]])])

    print('Running samtools:mpileup for %s' % basename, file=sys.stderr)
    out_vcf = os.path.join(output_folder, basename + '_mpileup_output')
    os.makedirs(out_vcf)

    for item in regions:
        chromosome, start, end, strand, bam_file, region_basename = item
        region = '%s%s%s%s%s' % (chromosome, ":", int(start), "-", int(end))
        output = os.path.join(out_vcf, region_basename + '.vcf')

        cl_vcf = 'samtools mpileup -v --region %s --fasta-ref %s %s > %s' % (region, reference_genome, bam_file, output)
        subprocess.check_call(cl_vcf, shell=True, env=os.environ.copy())

    print('Collecting variants for %s' % basename, file=sys.stderr)
    out_bcf = os.path.join(output_folder, basename + 'output_bcftools')
    os.makedirs(out_bcf)

    vcf_files = [f for f in os.listdir(out_vcf) if os.path.isfile(os.path.join(out_vcf, f))]
    for arch in vcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf'):
            name = arch[:-4]
            output = os.path.join(out_bcf, name + '_BCFcall.vcf')

            cl_bcf = 'bcftools call -v -c %s > %s' % (''.join([out_vcf, arch]), output)
            subprocess.check_call(cl_bcf, shell=True, env=os.environ.copy())

    print('Collecting significant variant calls for %s' % basename, file=sys.stderr)
    out_svc = os.path.join(output_folder, basename + 'output_svc')
    os.makedirs(out_svc)

    bcf_files = [f for f in os.listdir(out_bcf) if os.path.isfile(os.path.join(out_bcf, f))]
    for arch in bcf_files:
        if not arch.startswith('.') and arch.endswith('.vcf'):
            name = arch[:-12]
            output = os.path.join(out_bcf, name + '_SIGNFcall.txt')

            cl_sed = "sed -n '/##/!p' %s | awk 'FNR>1' > %s" % (''.join([out_bcf, arch]), output)
            subprocess.check_call(cl_sed, shell=True, env=os.environ.copy())

    print('Consolidating all the significant variant calls for %s' % basename, file=sys.stderr)
    header = ['targetsite', 'site_name', 'chromosome', 'one_based_position', 'reference', 'variant', 'quality', 'genotype', 'depth', 'PL']
    variants = list()

    svc_files = [f for f in os.listdir(out_svc) if os.path.isfile(os.path.join(out_svc, f))]
    for arch in svc_files:
        if not arch.startswith('.') and arch.endswith('.txt'):
            tag = arch[:-14]
            f = open(''.join([out_svc, arch]), 'r')
            reads = f.readlines()
            f.close()

            for line in reads:
                item = line.split()
                if 'INDEL' in item[7]:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[3][3:]] + ['_'.join(item[9][4:].split(','))])
                else:
                    variants.append(
                        [basename, tag] + item[:2] + item[3:6] + [str(int(item[9][0])) + '|' + str(int(item[9][2]))] +
                        [item[7].split(';')[0][3:]] + ['_'.join(item[9][4:].split(','))])

    out_file = open(out_snp, 'w')
    print(*header, sep='\t', file=out_file)
    for item in variants:
        print(*item, sep='\t', file=out_file)
    out_file.close()

    return variants


"""
Obtain variant off-target sequences
"""
def realignVariantBulge(bulge_sequence, window_sequence_variant, bulge_strand):
    bseq = bulge_sequence.replace('-', '')
    if bulge_strand == '+':
        m_bulge = re.search(bseq, window_sequence_variant, re.I)
    else:
        m_bulge = re.search(bseq, reverseComplement(window_sequence_variant), re.I)
    variant_bseq = m_bulge.group()
    variant_bseq = variant_bseq[:bulge_sequence.find('-')] + '-' + variant_bseq[bulge_sequence.find('-'):]

    return variant_bseq


def SNPreader(snp_file, min_quality=20, max_reference=10, min_depth=3):
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode='O')

    for snp in snp_file:
        basename, snpID, chromosome, one_based_position, reference, variant, quality, genotype, depth, PL = snp

        if (float(quality) >= int(min_quality)) and (len(reference) <= int(max_reference)) and (int(depth) >= min_depth):
            targetsite = snpID.split(':')[0]
            position = int(one_based_position) - 1
            key = '_'.join([targetsite, chromosome])
            ga[HTSeq.GenomicInterval(chromosome, position, position + 1, ".")] = [position, reference, variant, genotype, key]

    return ga


def arrayOffTargets(matched_file, search_radius=20):
    offtargets_dict = {}
    gi_dict = {}

    with open(matched_file, 'r') as g:
        g.readline()
        for line in g:
            site = line.strip().split('\t')

            Chromosome = site[0]
            start = int(site[1]) - search_radius
            end = int(site[2]) + search_radius
            Name = site[3]

            offtargets_dict[Name] = site

            gi_dict[Name] = HTSeq.GenomicInterval(Chromosome, start, end, ".")

    return offtargets_dict, gi_dict


def snpAdjustment(matched_file, snp_file, out, mismatch_threshold, min_quality=25, max_reference=5, min_depth=5, search_radius=20):
    output_file = open(out, 'w')
    print('Chromosome', 'PositionStart', 'PositionEnd', 'Name', 'ReadCount',
          'WindowName', 'WindowChromosome', 'Variant_WindowSequence',
          'Variant_Site_SubstitutionsOnly.Sequence', 'Variant_Site_SubstitutionsOnly.NumSubstitutions',
          'Variant_Site_SubstitutionsOnly.Strand', 'Variant_Site_SubstitutionsOnly.Start', 'Variant_Site_SubstitutionsOnly.End',
          'Variant_Site_GapsAllowed.Sequence', 'Variant_Site_GapsAllowed.Length', 'Variant_Site_GapsAllowed.Score',
          'Variant_Site_GapsAllowed.Substitutions', 'Variant_Site_GapsAllowed.Insertions', 'Variant_Site_GapsAllowed.Deletions',
          'Variant_Site_GapsAllowed.Strand', 'Variant_Site_GapsAllowed.Start', 'Variant_Site_GapsAllowed.End',
          'FileName', 'Cell', 'Targetsite', 'FullName', 'TargetSequence', 'Variant_RealignedTargetSequence',
          'Reference', 'Variant', 'Genotype',
          sep='\t', file=output_file)
    output_file.close()

    offtargets, gi_offtargets = arrayOffTargets(matched_file, search_radius)
    ga_snp = SNPreader(snp_file, min_quality, max_reference, min_depth)

    for name in offtargets:
        site = offtargets[name]
        gi = gi_offtargets[name]

        chromosome = site[0]
        PositionStart = int(site[1])
        PositionEnd = int(site[2])
        window_sequence = site[7]
        window_sequence = window_sequence.upper()

        fileName, cell, targetsite, fullName, TargetSequence = site[22:27]
        output01 = site[0:7]
        output03 = [fileName, cell, targetsite, fullName, TargetSequence]

        #  obtain variant window sequence
        wkey = '_'.join([targetsite, chromosome])

        insert_start, insert_end, insert_var, snp_data = list(), list(), list(), {}

        for i, v in ga_snp[gi].steps():
            if v:
                position, reference, variant, genotype, key = v
                if key == wkey:
                    variant = variant.split(',')[0]
                    for n, pos in enumerate(range(gi.start, gi.end)):
                        if pos == int(position):
                            insert_var.append(variant.lower())
                            insert_start.append(n)
                            end_pos = n + len(reference)
                            insert_end.append(end_pos)
                            snp_data[str(position)] = [position, reference, variant, genotype]

        tri = 0
        window_sequence_variant = ''
        for i in range(len(insert_var)):
            variant = insert_var[i]
            pos = insert_start[i]
            window_sequence_variant += window_sequence[tri:pos] + variant.lower()
            tri = insert_end[i]
        window_sequence_variant += window_sequence[tri:]

        #  get new off-target sequences
        window_sequence_var = window_sequence_variant.upper()

        #  only proceed if there is a variant in the window sequence
        if window_sequence_variant != window_sequence:
            offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, \
            chosen_alignment_strand_b, bulged_start, bulged_end, realigned_target = \
                alignSequences(TargetSequence, window_sequence_var, mismatch_threshold=mismatch_threshold)

            # get genomic coordinates of sequences
            mm_start, mm_end, b_start, b_end = '', '', '', ''
            if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '+':
                mm_start = PositionStart - search_radius + int(start_no_bulge)
                mm_end = PositionStart - search_radius + int(end_no_bulge)
            if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '-':
                mm_start = PositionEnd + search_radius - int(end_no_bulge)
                mm_end = PositionEnd + search_radius - int(start_no_bulge)

            if bulged_offtarget_sequence and chosen_alignment_strand_b == '+':
                b_start = PositionStart - search_radius + int(bulged_start)
                b_end = PositionStart - search_radius + int(bulged_end)
            if bulged_offtarget_sequence and chosen_alignment_strand_b == '-':
                b_start = PositionEnd + search_radius - int(bulged_end)
                b_end = PositionEnd + search_radius - int(bulged_start)

            # collect snp data from actual variant sequences
            total_genotype, total_reference, total_variant = '', '', ''
            bk_mm, bk_bulge = [], []
            variant_ots_no_bulge, variant_ots_bulge = '', ''
            for pos in snp_data:
                position, reference, variant, genotype = snp_data[pos]
                if offtarget_sequence_no_bulge:
                    if position >= mm_start and position < mm_end and position not in bk_mm:
                        bk_mm.append(position)
                        if total_genotype != '':
                            total_genotype += ''.join([':', genotype])
                            total_reference += ''.join([':', reference])
                            total_variant += ''.join([':', variant])
                        else:
                            total_genotype += ''.join([genotype])
                            total_reference += ''.join([reference])
                            total_variant += ''.join([variant])
                        if chosen_alignment_strand_m == '+':
                            m_no_bulge = re.search(offtarget_sequence_no_bulge, window_sequence_variant, re.I)
                        else:
                            m_no_bulge = re.search(offtarget_sequence_no_bulge, reverseComplement(window_sequence_variant), re.I)
                        variant_ots_no_bulge = m_no_bulge.group()

                if bulged_offtarget_sequence:
                    if position >= b_start and position < b_end and position not in bk_bulge:
                        bk_bulge.append(position)
                        if total_genotype != '':
                            total_genotype += ''.join([':', genotype])
                            total_reference += ''.join([':', reference])
                            total_variant += ''.join([':', variant])
                        else:
                            total_genotype += ''.join([genotype])
                            total_reference += ''.join([reference])
                            total_variant += ''.join([variant])
                        variant_ots_bulge = realignVariantBulge(bulged_offtarget_sequence, window_sequence_variant, chosen_alignment_strand_b)

                if not offtarget_sequence_no_bulge and not bulged_offtarget_sequence:
                    if total_genotype != '':
                        total_genotype += ''.join([':', genotype])
                        total_reference += ''.join([':', reference])
                        total_variant += ''.join([':', variant])
                    else:
                        total_genotype += ''.join([genotype])
                        total_reference += ''.join([reference])
                        total_variant += ''.join([variant])

            #  only continue if there are changes in the off-target sequences
            if total_genotype:
                output02 = [variant_ots_no_bulge, mismatches, chosen_alignment_strand_m, str(mm_start), str(mm_end),
                            variant_ots_bulge, length, score, substitutions, insertions, deletions,
                            chosen_alignment_strand_b, str(b_start), str(b_end)]
                output04 = [total_reference, total_variant, total_genotype]
                output_line = output01 + [window_sequence_variant] + output02 + output03 + [realigned_target] + output04
                with open(out, 'a') as output_file:
                    print(*output_line, sep='\t', file=output_file)


"""
Main function
"""
def callVariants(matched_file, ref, bam_file, out, search_radius, mismatch_threshold, min_quality=20, max_reference=10, min_depth=3):
    reference_genome = pyfaidx.Fasta(ref)
    basename = os.path.basename(matched_file)

    output_folder = os.path.dirname(out)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    out_snp = os.path.join(output_folder, basename, '_mpileupCall.txt')
    snp_file = snpCall(basename, matched_file, reference_genome, bam_file, out_snp, search_radius)

    snpAdjustment(matched_file, snp_file, out, mismatch_threshold, min_quality, max_reference, min_depth, search_radius)


def main():
    parser = argparse.ArgumentParser(description='Implement samtools:mpileup to identify genomic variants and adjust the off-target sequence when required.')
    parser.add_argument("--matched_file", help="full_path_to/matched file in 'identified' folder", required=True)
    parser.add_argument('--ref', help="Reference Genome Fasta", required=True)
    parser.add_argument('--bam', help="Sorted BAM file", required=True)
    parser.add_argument('--search_radius', help="Search radius around the position window", default=20, type=int)
    parser.add_argument('--mismatch_threshold', help='Maximum score threshold', default=7, type=int)
    parser.add_argument('--minimum_quality', help="Minimum quality score from samtools:mpileup threshold", default=20, type=int)
    parser.add_argument("--max_reference", help="Maximum reference size threshold", default=10)
    parser.add_argument("--min_depth", help="Minimum number of reads from which the SNP assessment was taken.", default=3)
    parser.add_argument('--out', help="full_path_to/output_file", required=True)
    args = parser.parse_args()

    callVariants(args.matched_file, args.ref, args.bam, args.out, args.search_radius, args.mismatch_threshold, args.minimum_quality, args.max_reference, args.min_depth)

if __name__ == "__main__":
    main()
