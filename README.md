# sequencingData_track_R
The R code was developed to draw the genomic track of signals for next-generation sequencing data. 

Usage: Rscript R.browser.r [options]

Options:
        -d CHARACTER, --directory=CHARACTER
                To specify a directory which stores BEDGRAPH files. Default: NULL. Required.

        -o CHARACTER, --out_directory=CHARACTER
                To specify an output directory. Default: the same with directory.

        -m CHARACTER, --mode=CHARACTER
                To specify a way for single task or a batch of tasks (single or batch). Default: single.

        -f CHARACTER, --filename=CHARACTER
                To specify a filenames to store a series of genomic coordinates, such as: 'chrom start end...' by TAB. Default: NULL.

        -c CHARACTER, --chrom=CHARACTER
                To specify a chromosome name. Default: NULL.

        -s CHARACTER, --start=CHARACTER
                To specify a start of genomic coordinate. Default: NULL.

        -e CHARACTER, --end=CHARACTER
                To specify an end of genomic coordinate. Default: NULL.

        -p CHARACTER, --pattern=CHARACTER
                To specify a string pattern embedding in file names of BEDGRAPH files. Default: bdg.

        -n CHARACTER, --name=CHARACTER
                To specify a string pattern in the output figures. Default: Region.chrom_start_end.

        --show_ylim=CHARACTER
                To specify whether to show ylim. Default: TRUE for each.

        --is_commmon_ylim=CHARACTER
                To specify whether to draw the same scale of ylim. Default: FALSE.

        --is_bigData=CHARACTER
                To specify whether to use AWK mode. Default: TRUE.

        --layout_mode=CHARACTER
                To specify which layout mode will be used (speparate or collapse). Default: speparate.

        --scale_show=CHARACTER
                To specify the length to show in visualization figure. Default: 5000 (that is 5kb).

        --cex=CHARACTER
                To specify the size to show labels. Default: 3.5.

        --figure_width=CHARACTER
                To specify the width of figure. Default: 2500.

        --figure_height=CHARACTER
                To specify the height of figure. Default: 600.

        --is_smooth=CHARACTER
                To specify whether to smooth read signals. Default: FALSE.

        --labels=CHARACTER
                To specify the labels to show; format: lab1,lab2,...labN. Default: the labels are splitted by 'PATTEN'.

        --colors=CHARACTER
                To specify the colors corresponding to labels; format: color1,color2,...colorN. Default: colors will be grouped by labels (split by 'rep').

        --species=CHARACTER
                To specify the species name to show. Default: mm10.

        --colorPattern=CHARACTER
                To specify a pattern to split file names to group them for showing. Default:'.*/|rep.*'.

        --color_transparent=CHARACTER
                To specify the degree of transparent color between 0 (transparent) and 1. Default: 1.

        --layout_width=CHARACTER
                To specify the ratio of the label panel over the figure panel. Default: '1:10'.

        --scale_signal=CHARACTER
                To specify an integer factor to scale the read signal data. Default: '1:10'.

        --layout_config=CHARACTER
                To specify the height for each track. Default: 3 for each.

        -h, --help
                Show this help message and exit
                
For example: 
1. Linux/macOS
Rscript $HOME/Visual/R.browser.r -d $HOME/Visual -m batch -f $HOME/Visual/Final.list.txt -p .bdg --labels P6_WT,P23_WT,P6_Math5KO,P23_Math5KO.DMR --colors blue,blue,red,red,black --show_ylim T,T,T,T,F --layout_config 3,3,3,3,1.5 --color_transparent 1 --figure_width 600 --figure_height 300 --layout_width 1:3 -o $HOME/Visual
2. Windows
/Progra~1/R/R-3.5.0/bin/Rscript.exe ./Figures/Visual/R.browser.r -d Figures/Visual -m batch -f ./Figures/Visual/Final.list.txt -p .bdg --labels P6_WT,P23_WT,P6_Math5KO,P23_Math5KO,DMR --colors blue,blue,red,red,black --show_ylim T,T,T,T,F --layout_config 3,3,3,3,1.5 --color_transparent 1 --figure_width 600 --figure_height 300 --layout_width 1:3 -o ./Figures/Visual
