# this script confirms that the between-species 
# orthologs produced by mapPeaks.py match Alan's
# list of orthologs (venn diagram)

# quick name repare function
def namefix(string):
	string = string.replace("Rpd3", "HDAC1")
	string = string.replace("Eip74EF", "Eip74")
	return string.upper()


# Species orthologs from Alan Boyle's venn diagram:
hs_boyleHsCe = map(namefix, ["POLR2A", "POLR3A", "MEF2A", "MEF2C", "PAX5", "ZEB1", "NR2F2", "GATA2", "GATA3", "FOXA1", "FOXA2", "SIX5", "E2F4", "RXRA", "E2F1", "NR4A1", "PRDM1", "E2F6", "NR2C2", "ESRRA", "MXI1", "FOSL1", "HNF4G", "FOS", "MYBL2", "FOSL2", "HNF4A", "ATF3"])
ce_boyleHsCe = map(namefix, ["AMA-1", "RPC-1", "MEF-2", "MEF-2", "PAX-1", "ZAG-1", "UNC-55", "ELT-3", "ELT-3", "PHA-4", "PHA-4", "UNC-39", "EFL-1", "UNC-55", "EFL-1", "NHR-6", "BLMP-1", "EFL-1", "NHR-67", "NHR-6", "MDL-1", "FOS-1", "NHR-25", "FOS-1", "GEI-11", "FOS-1", "NHR-25", "FOS-1"])


hs_boyleHsDm = map(namefix, ["HDAC2", "HDAC1", "POLR2A", "TBP", "HDAC2", "HDAC1", "SIN3A", "RXRA", "HDAC8", "HDAC8", "NR2F2", "STAT5A", "CTCF", "CTCFL", "TRIM28", "NR2C2", "NR2C2", "STAT3", "STAT1", "STAT2", "HDAC6", "ELF1"])
dm_boyleHsDm = map(namefix, ["Rpd3", "Rpd3", "RpII215", "Tbp", "Hdac3", "Hdac3", "Sin3A", "usp", "Rpd3", "Hdac3", "usp", "Stat92E", "CTCF", "CTCF", "bon", "tll", "Hr78", "Stat92E", "Stat92E", "Stat92E", "HDAC4a", "Eip74EF"])


ce_boyleCeDm = map(namefix, ["AMA-1", "LIN-39", "MAB-5", "UNC-62", "CEH-28", "UNC-55", "NHR-67", "NHR-67", "NHR-10", "NHR-2"])
dm_boyleCeDm = map(namefix, ["RpII215", "Dfd", "Dfd", "hth", "Dll", "usp", "tll", "Hr78", "kni", "kni"])



# Species orthologs produced with mapPeaks.py:
hs_orthoHsCe = map(namefix, ["PRDM1", "MYBL2", "E2F1", "E2F4", "E2F6", "PAX5", "FOXA2", "FOXA1", "POLR3A", "FOSL1", "FOSL2", "ATF3", "FOS", "MEF2A", "MEF2C", "NR4A1", "ESRRA", "HNF4G", "HNF4A", "MXI1", "SIX5", "GATA3", "GATA2", "ZEB1", "NR2C2", "NR2F2", "RXRA", "POLR2A"])
ce_orthoHsCe = map(namefix, ["BLMP-1", "GEI-11", "EFL-1", "PAX-1", "PHA-4", "RPC-1", "FOS-1", "MEF-2", "NHR-6", "NHR-25", "MDL-1", "UNC-39", "ELT-3", "ZAG-1", "NHR-67", "UNC-55", "AMA-1"])


hs_orthoHsDm = map(namefix, ["STAT1", "STAT2", "STAT3", "STAT5A", "HDAC6", "CTCFL", "CTCF", "ELF1", "TBP", "TRIM28", "HDAC8", "HDAC2", "HDAC1", "SIN3A", "POLR2A", "NR2C2", "NR2F2", "RXRA"])
dm_orthoHsDm = map(namefix, ["Stat92E", "HDAC4a", "CTCF", "Eip74", "TBP", "bon", "HDAC3", "HDAC1", "Sin3A", "RpII215", "tll", "Hr78", "USP"])


ce_orthoCeDm = map(namefix, ["LIN-39", "MAB-5", "UNC-62", "NHR-10", "NHR-2", "CEH-28", "AMA-1", "NHR-67", "UNC-55"])
dm_orthoCeDm = map(namefix, ["Dfd", "Hth", "KNI", "Dll", "RpII215", "tll", "Hr78", "USP"])


print
print "HS/CE:"
b1, b2 = hs_boyleHsCe, ce_boyleHsCe
o1, o2 = hs_orthoHsCe, ce_orthoHsCe
diff1bo = set(b1).difference(set(o1))
diff1ob = set(b1).difference(set(o1))
diff2bo = set(b2).difference(set(o2))
diff2ob = set(o2).difference(set(b2))
print diff1bo, diff1ob
print diff2bo, diff2ob
print len(set(b1)), len(set(o1))
print len(set(b2)), len(set(o2))

print
print "HS/DM:" 
b1, b2 = hs_boyleHsDm, dm_boyleHsDm
o1, o2 = hs_orthoHsDm, dm_orthoHsDm
diff1bo = set(b1).difference(set(o1))
diff1ob = set(b1).difference(set(o1))
diff2bo = set(b2).difference(set(o2))
diff2ob = set(o2).difference(set(b2))
print diff1bo, diff1ob
print diff2bo, diff2ob
print len(set(b1)), len(set(o1))
print len(set(b2)), len(set(o2))

print
print "CE/DM:" 
b1, b2 = ce_boyleCeDm, dm_boyleCeDm
o1, o2 = ce_orthoCeDm, dm_orthoCeDm
diff1bo = set(b1).difference(set(o1))
diff1ob = set(b1).difference(set(o1))
diff2bo = set(b2).difference(set(o2))
diff2ob = set(o2).difference(set(b2))
print diff1bo, diff1ob
print diff2bo, diff2ob
print len(set(b1)), len(set(o1))
print len(set(b2)), len(set(o2))


#python ~/Desktop/orthologComparison.py