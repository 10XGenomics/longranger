# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
import tenkit.bio_io as tk_io

class Event(object):
    def __init__(self):
        #self.set_filters()
        self.name ="."
        self.filter_rst = (True, "PASS")

    def set_filters(self, Flt_pval = 1e-19, Flt_WildCov = 6.4, Flt_BGImparityPval = 3e-3, Flt_BGBaseCov = 6.2, Flt_BGPhaseFrac = 0.52):
        self.Flt_pval = Flt_pval
        self.Flt_WildCov = Flt_WildCov
        self.Flt_BGImparityPval = Flt_BGImparityPval
        self.Flt_BGBaseCov = Flt_BGBaseCov
        self.Flt_BGPhaseFrac = Flt_BGPhaseFrac
        self.filter_rst = self.isPassFilter()

    @staticmethod
    def from_native(line):
        ev = Event()
        ev.read_from_native(line)
        return ev

    # Load from 'native' format created by call_cnv
    def read_from_native(self, line, ZS="HET"):
        self.ZS = ZS
        info = line.split()

        self.hasBGSanityMetrics = False
        self.hasDetails = False

        if len(info) >= 21:
            self.hasBGSanityMetrics = True
        if len(info) > 4:
            self.hasDetails = True

        self.chrom=info[0]
        self.start=int(info[1])
        self.end=int(info[2])
        if self.start > self.end: 
            self.start, self.end = self.end, self.start

        if self.hasDetails:
            self.pval=float(info[4])
            self.hap = int(info[5])
            self.WildCov = int(info[6])
            self.size = int(info[7])
            self.MolTtl = int(info[8])
            self.MolTtlNoR = int(info[9])
            self.MolWild = int(info[10])
            self.MolWildNoR = int(info[11])
            self.MolDel = int(info[12])
            self.MolDelNoR = int(info[13])

        if self.hasBGSanityMetrics:
            self.BGSize = int(info[14])
            self.BGImparityPval = float(info[15])
            self.BGTtlRCnt = int(info[16])
            self.BGHP1RCnt = int(info[17])
            self.BGHP2RCnt = int(info[18])
            self.BGBaseCov = float(info[19])
            self.BGPhaseFrac = float(info[20])

    @staticmethod
    def from_bedpe(line):
        ev = Event()
        ev.read_from_bedpe(line)
        return ev

    # Load from BEDPE
    def read_from_bedpe(self, line):
        fields = line.split()

        # FIXME -- read HAP field

        self.chrom=fields[0]
        self.start=int(fields[1])
        self.end=int(fields[4])
        #self.pval=float(fields[7])
        self.filters = fields[8]
        info_str = fields[11]

        info_pairs = info_str.strip(";").split(";")
        info = {}
        for pair in info_pairs:
            parts = pair.split("=");
            info[parts[0]] = parts[1]

        self.ZS = info["ZS"]
        self.HAP = info["HAP"]
        self.WildCov = float(info["WildCov"])
        self.size = int(info["size"])
        self.MolTtl = int(info["MolTtl"])
        self.MolTtlNoR = int(info["MolTtlNoR"])
        self.MolDel = int(info["MolDel"])
        self.MolDelNoR = int(info["MolDelNoR"])
        self.MolWild = int(info["MolWild"])
        self.MolWildNoR = int(info["MolWildNoR"])
        self.pval = float(info["PVAL"])

        self.BGSize = int(info["BGSize"])
        self.BGImparityPval = float(info["BGImparityPval"])
        self.BGTtlRCnt = int(info["BGTtlRCnt"])
        self.BGHP1RCnt = int(info["BGHP1RCnt"])
        self.BGHP2RCnt = int(info["BGHP2RCnt"])
        self.BGBaseCov = float(info["BGBaseCov"])
        self.BGPhaseFrac = float(info["BGPhaseFrac"])

    # Apply filters to event
    def isPassFilter(self):
        if not self.hasDetails:
            return (True, "PASS")
        elif self.pval > self.Flt_pval:
            return (False, "LARGE_PVALUE")
        elif self.WildCov < self.Flt_WildCov:
            return (False, "LOW_WILD_COV")
        elif self.BGSize ==0:
            return (True, "PASS")
        elif self.BGImparityPval < self.Flt_BGImparityPval:
            return (False, "IMPARITY_IN_HAP_COV")
        elif self.BGBaseCov < self.Flt_BGBaseCov:
            return (False, "LOW_COV_REGION")
        elif self.BGPhaseFrac < self.Flt_BGPhaseFrac:
            return (False, "LOW_PHASING_REGION")
        else:
            return (True, "PASS")
        
        
    def check_blacklist(self, blacklist, blacklist_name):
        if blacklist and self.chrom in blacklist and blacklist[self.chrom].overlapping_regions(self.start, self.end):
            if self.filter_rst[0]: 
                self.filter_rst= (False, blacklist_name)
            else: 
                self.filter_rst = (False, self.filter_rst[1]+";"+blacklist_name)

    def check_whitelist(self, whitelist, whitelist_name):
        if not self.chrom in whitelist or \
            not whitelist[self.chrom].overlapping_regions(self.start, self.end):
            if self.filter_rst[0]:
                self.filter_rst= (False, whitelist_name)
            else:
                self.filter_rst = (False, self.filter_rst[1]+";"+whitelist_name)

    # Legacy output format
    def get_bed_output(self):
        passString = self.filter_rst[1]
        region = "%s\t%d\t%d" % (self.chrom, self.start, self.end)
        OutString =region+"\t"+passString

        if self.hasBGSanityMetrics:
            basic = "\t%.4g;%.4g;%d;%d;%d;%d;%d;%d;%d" % \
            (self.pval, self.WildCov, self.size, self.MolTtl, self.MolTtlNoR, self.MolWild, self.MolWildNoR, self.MolDel, self.MolDelNoR)
            OutString += basic

        if self.hasBGSanityMetrics:
            check = "\t%d;%.4g;%d;%d;%d;%.4g;%.4g" % \
            (self.BGSize, self.BGImparityPval, self.BGTtlRCnt, self.BGHP1RCnt, self.BGHP2RCnt, self.BGBaseCov, self.BGPhaseFrac)
            OutString += check

        return OutString

    # BEDPE output
    def get_bedpe_output(self):
        part1 = "%s\t%d\t%d\t%s\t%d\t%d\t" % (self.chrom, self.start, self.start+1, self.chrom, self.end, self.end+1)

        filters = self.filter_rst[1]

        # INFO fields
        info = "ZS=%s;HAPS=%d,%d;SOURCE=CNV;" % (self.ZS, self.hap, self.hap)

        info += "WildCov=%.4g;Size=%d;MolTtl=%d;MolTtlNoR=%d;MolDel=%d;MolDelNoR=%d;MolWild=%d;MolWildNoR=%d;PVAL=%.2g;" % \
                (self.WildCov, self.size, self.MolTtl, self.MolTtlNoR, self.MolDel, self.MolDelNoR, self.MolWild, self.MolWildNoR, self.pval)

        info += "BGSize=%d;BGImparityPval=%.4g;BGTtlRCnt=%d;BGHP1RCnt=%d;BGHP2RCnt=%d;BGBaseCov=%.4g;BGPhaseFrac=%.4g;TYPE=DEL" % \
                (self.BGSize, self.BGImparityPval, self.BGTtlRCnt, self.BGHP1RCnt, self.BGHP2RCnt, self.BGBaseCov, self.BGPhaseFrac)

        #qual = int(min(1000, -10 * np.log10(max(self.pval, 1e-199))))
        #qual = self.MolDelNoR
        part2 = "%s\t%d\t.\t.\t%s\t%s" % (self.name, self.MolDelNoR, filters, info)
        return part1 + part2

    
def do_blacklist_filtering(events, blacklist_map):
    for bl_name, bl_list in blacklist_map.iteritems():
        with open(bl_list) as fBL:
            blacklist = tk_io.get_target_regions(fBL)
            
            for e in events:
                e.check_blacklist(blacklist, bl_name)
    
def do_whiltelist_filtering(events, whitelist_file, whitelist_name):
    with open(whitelist_file) as fWT:
        whitelist = tk_io.get_target_regions(fWT)
        
        for e in events:
            e.check_whitelist(whitelist, whitelist_name)
    

class HomoEvent(object):
    def __init__(self):
        self.name = "."

    @staticmethod
    def from_native(line):
        ev = HomoEvent()
        ev.read_from_native(line)
        return ev


    def read_from_native(self, line):
        fields = line.split()

        self.chrom=fields[0]
        self.start=int(fields[1])
        self.end=int(fields[2])
        self.filters = fields[3]
        self.pval = float(fields[4])
        self.Nread = float(fields[5])
        self.ttlMolMs = int(fields[6])
        if self.start > self.end: 
            self.start, self.end = self.end, self.start


    @staticmethod
    def from_bedpe(line):
        ev = HomoEvent()
        ev.read_from_bedpe(line)
        return ev


    def read_from_bedpe(self,line):
        fields = line.split()

        self.chrom=fields[0]
        self.start=int(fields[1])
        self.end=int(fields[4])
        self.name=fields[6]
        self.ttlMolMs = int(fields[7])
        self.filters = fields[8]

        info_str = fields[11]
        info_pairs = info_str.strip(";").split(";")
        info = {}
        for pair in info_pairs:
            parts = pair.split("=");
            info[parts[0]] = parts[1]

        self.Nread = info["NRead"]
        self.ZS = info["ZS"]
        self.HAPS = info["HAPS"]
        self.pval = info["PVAL"]
        
    def check_blacklist(self, blacklist, blacklist_name):
        if blacklist and self.chrom in blacklist and \
            blacklist[self.chrom].overlapping_regions(self.start, self.end):
            self.filters = blacklist_name
        else:
            self.filters = "PASS"

    def check_whitelist(self, whitelist, whitelist_name):
        if whitelist and self.chrom in whitelist and \
            whitelist[self.chrom].overlapping_regions(self.start, self.end):
            self.filters = "PASS"
        else:
            self.filters = whitelist_name

    
    def get_bedpe_output(self):
        part1 = "%s\t%d\t%d\t%s\t%d\t%d\t" % (self.chrom, self.start, self.start+1, self.chrom, self.end, self.end+1)

        # INFO fields
        info = "ZS=HOM;HAPS=1,1;SOURCE=CNV;Size=%d;NRead=%d;PVAL=%.2g;" % (self.end - self.start, self.Nread, self.pval)

        #qual = int(min(1000, -10 * np.log10(self.pval)))
        if self.ttlMolMs <= 0:
            part2 = "%s\t%d\t.\t.\t%s\t%s" % (self.name, 10, self.filters, info)
        else:
            part2 = "%s\t%d\t.\t.\t%s\t%s" % (self.name, self.ttlMolMs, self.filters, info)
        return part1 + part2

    def get_bed_output(self):
        return "%s\t%d\t%d\t%s\t%.4g\t%.1f\t%d" % (self.chrom, self.start, self.end, self.filters, self.pval, self.Nread, self.ttlMolMs)


def do_homo_whiltelist_filtering(events, whitelist_file, whitelist_name):
    with open(whitelist_file) as fWT:
        whitelist = tk_io.get_target_regions(fWT)
        
        for e in events:
            e.check_whitelist(whitelist, whitelist_name)
    



class SimpleEvent(object):
        def __init__(self, line):
                fields = line.split()

                self.chrom=fields[0]
                self.start=int(fields[1])
                self.end=int(fields[2])
                self.isPass = fields[3]=="PASS"
