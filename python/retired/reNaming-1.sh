#!/usr/bin/sh



# Renaming scheme (Part 1):  hs_selection > hs_standards (local)

# Rename IDR peak sets:
#reNamr ~/meTRN/idr/peaks/ hs_selection hs_standards

# Rename peak sets:
#reNamr ~/meTRN/data/peaks/ hs_selection hs_standards
#reNamr ~/meTRN/data/peaks/ orthoHs standHs

# Rename HOT simulations:
#reNamr ~/meTRN/data/hot/simulation/montecarlo/ hs_selection hs_standards
#reNamr ~/meTRN/data/hot/simulation/analysis/ hs_selection hs_standards

# Rename HOT regions:
#reNamr ~/meTRN/data/hot/analysis/ hs_selection hs_standards
#reNamr ~/meTRN/data/hot/compare/ hs_selection hs_standards
#reNamr ~/meTRN/data/hot/density/ hs_selection hs_standards
#reNamr ~/meTRN/data/hot/regions/ hs_selection hs_standards
#reNamr ~/meTRN/data/hot/overlap/ hs_selection hs_standards

# Rename GO analysis:
#reNamr ~/meTRN/data/go/ hs_selection hs_standards
#reNamr ~/meTRN/data/go/hs_standards_com_cx_hct/p5e-1/results/ hs_selection hs_standards
#reNamr ~/meTRN/data/go/hs_standards_com_cx_hct/p5e-1/summary/ hs_selection hs_standards

# Rename MEME results:
#reNamr ~/meTRN/meme/ hs_selection hs_standards

# Rename Co-associations:
#reNamr ~/meTRN/data/coassociations/ orthoHs standHs
#reNamr ~/meTRN/data/coassociations/runs/ orthHs standHs
#reNamr ~/meTRN/data/coassociations/comparison/hs/ce/ orthoHs standHs

# Rename R graphs:
#reNamr ~/Desktop/Dropbox/meTRN/Code/R/graphs/ hs_selection hs_standards
#reNamr ~/Desktop/Dropbox/meTRN/Code/R/graphs/ orthoHs standHs



# Renaming scheme (Part 2):  hs_selection > hs_standards (remote)

# Rename peak sets:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/peaks/ hs_selection hs_standards
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/peaks/ orthoHs standHs

# Rename HOT simulations:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/hot/simulation/montecarlo/ hs_selection hs_standards
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/hot/simulation/analysis/ hs_selection hs_standards

# Rename HOT regions:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/hot/density/ hs_selection hs_standards

# Rename Co-associations:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/ orthoHs standHs
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/runs/ orthoHs standHs




# Renaming scheme (Part 3):  hs_alternate > hs_selection (local)

# Rename IDR peak sets:
#reNamr ~/meTRN/idr/peaks/ hs_alternate hs_selection

# Rename peak sets:
#reNamr ~/meTRN/data/peaks/ hs_alternate hs_selection

# Rename HOT simulations:
#reNamr ~/meTRN/data/hot/simulation/montecarlo/ hs_alternate hs_selection
#reNamr ~/meTRN/data/hot/simulation/analysis/ hs_selection hs_selection



# Renaming scheme (Part 4):  hs_alternate > hs_selection (remote)

# Rename peak sets:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/peaks/ hs_alternate hs_selection

# Rename HOT simulations:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/hot/simulation/montecarlo/ hs_alternate hs_selection
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/hot/simulation/analysis/ hs_alternate hs_selection

# Rename HOT regions:
#reNamr /srv/gs1/projects/snyder/claraya/meTRN/data/hot/density/ hs_alternate hs_selection
