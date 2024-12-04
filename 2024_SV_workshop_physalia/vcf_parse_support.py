#!/usr/bin/env python

import sys
import pysam


def main():
    if len(sys.argv) == 2:
        vcf = pysam.VariantFile(sys.argv[1])
    else:
        print("error VCF missing")
        sys.exit(0)
    for vcf_entry in vcf.fetch():
        suppv = vcf_entry.info.get("SUPP_VEC")
        print (suppv)


if "__main__" == __name__:
    main()
