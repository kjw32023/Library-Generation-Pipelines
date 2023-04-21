# load libraries
library(tidyverse)
library(Rsamtools)
library(seqinr)

# read in the table of edits
edits <- read_csv('MMR-responsive_+5GtoT_edit_candidates_v2_final.csv')

# read in IKW002 alignments
bam <- scanBam("lKW002_IKW002_NEB_hifi_stellar_cell_plasmid_pool_align.bam")

# names of the BAM fields
names(bam[[1]])
# [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"
# [9] "mrnm"   "mpos"   "isize"  "seq"    "qual"

# distribution of BAM flags
table(bam[[1]]$flag)

# function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)){
        structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
        do.call(c, x)
    }
}

# store names of BAM fields
bam_field <- names(bam[[1]])

# go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

# store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

# convert data frame to tibble for better compatibility with tidyverse functions
bam_df <- as_tibble(bam_df)

# remove the reverse read alignments
bam_df <- bam_df %>% filter(strand == '-')

logical <- logical(length=nrow(bam_df))

# add blank columns to bam_df
bam_df <- bam_df %>% mutate(guide_match = c(),
                            RT_match = c(),
                            target_match = c(),
                            no_recom = logical,
                            guide.RT_recom = logical,
                            RT.target_recom = logical,
                            guide.target_recom = logical)

# make an empty tibble
bam_rejoin <- tibble()

# the bowtie2 software align reads to pegRNA IDs by target sequence
## nrow(bam_rejoin) == sum(bam_rejoin$target_match)

# thus, the data assume the target site is never recombined
## target_match always == TRUE

# to ID data where recombination happened between the RT and target site,
# can't search for guide_match == TRUE & RT_match == TRUE
# because this is only possible if the correct target site is also there

# therefore, in this case only, recombination is identified by
# guide_match == FALSE & RT_match == FALSE

# find sequences that match the guide, PBS/RT, and target site
for (i in 1:nrow(edits)) {
    # get the guide, PBS/RT, target seqx, pegID
    pegID <- edits$pegRNA_ID[i]
    guide <-edits$Guide[i]
    RT <- edits$RT_template[i]
    target <- toupper(paste(rev(comp(unlist(strsplit(edits$Wide_Target_4A[i], ''))))))

    # search only strands that match the pegRNA IDs
    bam_df_peg <- bam_df %>% filter(rname == pegID)
    
    # search to see if the guide, RT, and target match the sequence
    bam_df_peg <- bam_df_peg %>% 
        mutate(guide_match = grepl(guide, bam_df_peg$seq, fixed=TRUE),
               RT_match = grepl(RT, bam_df_peg$seq, fixed=TRUE),
               target_match = grepl(target, bam_df_peg$seq, fixed=TRUE))

    for (j in 1:nrow(bam_df_peg)) {
        # no_recom: all() finds sites with no recombination
        bam_df_peg$no_recom[j] <- all(bam_df_peg$guide_match[j] == TRUE,
                                      bam_df_peg$RT_match[j] == TRUE,
                                      bam_df_peg$target_match[j] == TRUE)
        # rows with no recombination will satisfy all of the following code,
        # so eliminate them from further consideration
        if (bam_df_peg$no_recom[j] == FALSE) {
            # all() over RT, target finds recombination between guide/RT
            bam_df_peg$guide.RT_recom[j] <- all(bam_df_peg$guide_match[j] == FALSE,
                                                bam_df_peg$RT_match[j] == TRUE,
                                                bam_df_peg$target_match[j] == TRUE)
            # all() over guide, RT finds recombination between RT/target
            bam_df_peg$RT.target_recom[j] <- all(bam_df_peg$guide_match[j] == FALSE,
                                                 bam_df_peg$RT_match[j] == FALSE,
                                                 bam_df_peg$target_match[j] == TRUE)
            # all() over guide, target finds double recombination (between guide/target)
            bam_df_peg$guide.target_recom[j] <- all(bam_df_peg$guide_match[j] == TRUE,
                                                    bam_df_peg$RT_match[j] == FALSE,
                                                    bam_df_peg$target_match[j] == TRUE)
        }
    }
    # add pegRNAID specific table to new bam tibble
    bam_rejoin <- rbind(bam_rejoin, bam_df_peg)
}

# total # of aligned reads
tot <- nrow(bam_rejoin)
tot

# total # of aligned target sites (should == total # of reads)
tar_match <- sum(bam_rejoin$target_match)
tar_match

# total # of reads where guide matched alignment index
gui_match <- sum(bam_rejoin$guide_match)
gui_match

# total # of reads where RT matched alignment index
RT_match <- sum(bam_rejoin$RT_match)
RT_match

# total # of reads with mutations in both guide and RT
mut <- nrow(bam_rejoin) - (sum(bam_rejoin$guide_match) + sum(bam_rejoin$RT_match))
mut

# total # of reads with no recombination
no_recom <- sum(bam_rejoin$no_recom)
no_recom

# total # of reads recombined between guide and RT
gui.RT_recom <- sum(bam_rejoin$guide.RT_recom)
gui.RT_recom

# total # of reads recombined between RT and target (remove mutations)
RT.tar_recom <- sum(bam_rejoin$RT.target_recom) - mut
RT.tar_recom

# total # of doubly recombined reads (replacing RT between guide and target)
gui.tar_recom <- sum(bam_rejoin$guide.target_recom)
gui.tar_recom

# make a summary table
recombination <- c('No recombination', 'Guide/RT recombination',
                   'RT/Target recombination', 'Double recombination')
count <- c(no_recom, gui.RT_recom, RT.tar_recom, gui.tar_recom)
percent <- count/tot

summary <- tibble(recombination, count, percent)

# set column order
summary$recombination <- as.character(summary$recombination)
summary$recombination <- factor(summary$recombination, levels=unique(summary$recombination))

ggplot(summary, aes(x=recombination, y=percent)) +
    geom_col(fill='lightblue') +
    labs(x="Recombination Type",
         y ="Percent of
Alignments") +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color='lightgray', size=0.5, linetype = 1),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('recombination.svg', plot = last_plot(), device = 'svg', 
       path = 'Summary Figures', scale = 0.75, width = 6.5, height = 3, 
       units = "in", dpi = 300)
