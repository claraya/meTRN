#!/usr/bin/sh

# Backup the complete project structure:
#rsync -avz ~/meTRN /Volumes/Drobo/meTRN/


# Import ortholog families:
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/modENCODE/data/orthologs/family_group/* ~/meTRN/data/orthologs/families/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/modENCODE/data/orthologs/family_group/formatted* ~/meTRN/data/orthologs/families/
#reNamr ~/meTRN/data/orthologs/families/ cel ce
#reNamr ~/meTRN/data/orthologs/families/ dmel dm


# Transfer inputs, code, and qsub configurations:
#rsync -avz ~/meTRN/input/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/input/
#rsync -avz ~/meTRN/python/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/python/
#rsync -avz ~/meTRN/scripts/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/scripts/
#rsync -avz ~/meTRN/qsub/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/qsub/
#rsync -avz ~/Desktop/Dropbox/Code/Python/ claraya@crick.stanford.edu:~/Tools/python/


# Transfer peaks and HOT-region data:
#rsync -avz ~/meTRN/data/peaks/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/peaks/
#rsync -avz ~/meTRN/data/annotations/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/annotations/
#rsync -avz ~/meTRN/data/hot/density/mergeBed* claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/hot/density/
#rsync -avz ~/meTRN/data/hot/regions/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/hot/regions/


# Transfer peaks to modENCODE shared folder:
#rsync -avz ~/meTRN/data/peaks/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/modENCODE/data/peaks/sets/


# Transfer co-association results:
#rsync -avz araya@nexus.gs.washington.edu:~/meTRN/data/coassociations/*ortho* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/ce_selection* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/hs_selection* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/dm_selection* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/ce_corrected* ~/meTRN/data/coassociations/


# Transfer cellular resolution data (to Stanford SCG3):
#rsync -avz ~/meTRN/data/cells/peaks/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/cells/peaks/
#rsync -avz ~/meTRN/data/cells/annotations/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/cells/annotations/


# Transfer cellular resolution data (to UW Nexus):
#rsync -avz ~/meTRN/data/cells/peaks/ araya@nexus.gs.washington.edu:~/meTRN/data/cells/peaks/
#rsync -avz ~/meTRN/data/cells/annotations/ araya@nexus.gs.washington.edu:~/meTRN/data/cells/annotations/


# Transfer worm-specific co-association analyses:
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/ce_xptracked_reg_ex_all* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/ce_xptracked_reg_ex_240* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/ce_brunetlab* ~/meTRN/data/coassociations/


# Transfer worm-specific SOM data (Stanford):
#rsync -avz ~/meTRN/data/neurons/ claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/neurons/


# Recover worm-specific SOM seeds (Stanford):
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/neurons/ ~/meTRN/data/neurons/ 


# Pending:
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/*selection* ~/meTRN/data/coassociations/
#rsync -avz claraya@crick.stanford.edu:/srv/gs1/projects/snyder/claraya/meTRN/data/coassociations/*selection* ~/meTRN/data/coassociations/


# Transfer setup to GS clusters:
#rsync -avz ~/meTRN/input/ araya@nexus.gs.washington.edu:~/meTRN/input/
#rsync -avz ~/meTRN/python/ araya@nexus.gs.washington.edu:~/meTRN/python/
#rsync -avz ~/meTRN/scripts/ araya@nexus.gs.washington.edu:~/meTRN/scripts/
#rsync -avz ~/meTRN/qsub/ araya@nexus.gs.washington.edu:~/meTRN/qsub/
#rsync -avz ~/Tools/python/ araya@nexus.gs.washington.edu:~/Tools/python/

#rsync -avz ~/meTRN/data/peaks/*ortho* araya@nexus.gs.washington.edu:~/meTRN/data/peaks/
#rsync -avz ~/meTRN/data/annotations/ araya@nexus.gs.washington.edu:~/meTRN/data/annotations/
#rsync -avz ~/meTRN/data/neurons/ce_selection_com_ex_raw/binary/matrix/ araya@nexus.gs.washington.edu:~/meTRN/data/neurons/ce_selection_com_ex_raw/binary/matrix/
#rsync -avz ~/meTRN/data/neurons/ce_selection_com_cx_raw/binary/matrix/ araya@nexus.gs.washington.edu:~/meTRN/data/neurons/ce_selection_com_cx_raw/binary/matrix/


# Recover worm-specific SOM seeds (GS):
#rsync -avz araya@nexus.gs.washington.edu:~/meTRN/data/neurons/ce_selection_com_cx_raw/binary/reports/sam* ~/meTRN/data/neurons/ce_selection_com_cx_raw/binary/reports/



# Copy GO results for Alan:
#mkdir ~/Desktop/Dropbox/Shared\ with\ Alan/GO-tables
#cp ~/meTRN/data/go/ce_HOTregion_occ_cx_p05/p5e-1/summary/mapgo_complete_ce_HOTregion_occ_cx_p05_p5_hc1_hp5e-02_matrix ~/Desktop/Dropbox/Shared\ with\ Alan/GO-tables/
#cp ~/meTRN/data/go/dm_HOTregion_occ_cx_p05/p5e-1/summary/mapgo_complete_dm_HOTregion_occ_cx_p05_p5_hc1_hp5e-02_matrix ~/Desktop/Dropbox/Shared\ with\ Alan/GO-tables/
#cp ~/meTRN/data/go/hs_HOTregion_occ_cx_p05/p5e-1/summary/mapgo_complete_hs_HOTregion_occ_cx_p05_p5_hc1_hp5e-02_matrix ~/Desktop/Dropbox/Shared\ with\ Alan/GO-tables/


#top
