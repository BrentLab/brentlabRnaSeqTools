new                                # batchscript keyword new (new snapshot)
snapshotDirectory IGV_Snapshots    # output directory
maxPanelHeight 500                 # maximum height of igv browser viewer
genome organism.genome             # igv .genome file
load control.bam                   # load alignment files in bam_file_list
load perturbed.bam
goto chr1:1-100                    # load region of interest
snapshot batchfilename_locus1.png  # saves a snapshot of the IGV window to an image file
goto chr10:3-300                   # repeat at another locus
snapshot batchfilename_locus2.png
exit                               # quit session
