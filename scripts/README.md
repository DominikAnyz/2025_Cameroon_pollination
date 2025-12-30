
# Script organization

This README serves as a walkthrough of all the scripts contains in this
folder.

All the individual sripts have their own descriptions about what each
segment of code does, this should help with making sense of why
individual steps needed to be taken.

All analyses were carried out on MetaCentrum, which provides an online
computation infrastructure for researchers. Because of this, I was able
to run many of the bayesian models through rstan easily, however if you
are running the code on your computer, it might not work and certain
parameters for the individual models might have to be set differently.

COME BACK TO THIS!!!!!!!!!

## 01.Setup_for_analysis

The script contains all the packages, which are used throughout the
scripts. This is handy to run if you will be carrying out all the
analyses again, however is not necessary. If you are interested only in
a certain script, it might be easier to load the packages in each script
separately. If you do not have cmdstanr on your computer for running
Bayesian analyses, this might take a while.

## 02.Cleaning_raw_seedset_data

The script take the raw data collected from finished seeds and cleans
the dataset, so it dos not contain any unusable observations. The
dataset is striped of any seedsets, where the comment would suggest
infestation, loss of flower or uncertainty of treatment.

The document *Mt.Cameroon seed count raw.xlsx* is the data which was
collected. The only changes made to this document are that if a specimen
was missing a plant number, a unique number would be appointed to it as
to not be tied with any one other plant number.

The manipulation with the rest of the data was carried out using this
code in R

## 03.GLMM_pollination_indices

The script contains all the code used to run glmm’s for the four
hand-pollination indices, which were used to assess the effect of
elevation of them (natural seed set, pollen limitation, autogamy and
geitonogamy).

The script contains paths to saved models in order not to have to rerun
them again.

## 04.Setup_for_visitation_indices

The script is used to setup the datasets for the analyses of individual
pollination indices. It has a similar setup as
*06.Setup_for_index_comparisons*. This script has data aggregated at
individual observation level, whereas *06.Setup_for_index_comparisons*
has data aggregated at elevation and species level. Although it might
seem like a small change, the code would have been more complicated if
one script was to be made for both datasets, so we chose to separate it
into two.

Both of these setup scripts pull datasets *visitors.txt* and
*functional.txt*, which have been altered outside R. The process in
which the datasets were created is as follows:

The raw visitation data is in the form of excels, 7 for each species for
each elevation (provided that the species was in a given elevation).
These excels an be found in [📁
elevation_excels](/data/visitors/elevation_excels/).

It was necessary to combine all the excels into one. The steps for the
combining all the excels can be found in the document *Processing
steps.docx* in the folder **../data/visitors/process**.

TLDR; first it is necessary to check that all excels are filled in
correctly. This was done manually. All the excels can be found in the
folder **../data/visitors/elevation excels**. Next, it was necessary to
export every excel as a csv to make the files as small as possible. This
was done with the *XLS to CSV converter.hta* in the
**../data/visitors/process** folder. The csv’s (found in folder
**../data/visitors/csv files**) were then combined and the original
version of the combined document if under the name \*\_combined.csv\* in
the file **../data/visitors/csv files**.

### visitors.txt

This csv file was imported into an xls file labeled *Pollination
necessary.xlsx*. In this excel the unecessary columns were deleted.
Also, all the names of identified visitors were exported and added to
special excel *Cameroon grasslands species - morphospecies.xlsx* and it
was checked if any of the names used were not used differently between
the excels. The names were unified by the visitor identifiers (the
people who ID’d the visitors) every time one species was the same as
another, it is recorded in the excel and then the name for each given
species is bolded. Identifiers also recorded whether or not the visitor
can be identified as a “morphospecies” or not.

The information from this excel was transported into the excel
*Pollination necessary.xlsx*, from which the file *visitors.txt* was
created in folder **../data**.

### functional.txt

The text document *functional.txt* was created from the excel *Cameroon
grassland species - morphospecies.xlsx*. After all the visitors have
been identified, only the columns “SD.s.ID” (Sylvain Delabye’s ID) and
“functional.group” (which functional group the visitor belongs to) were
exported into the file *functional.txt* in folder **../data**

The following functional groups were used: Bee, Wasp, Hoverfly, Other
fly, Beetle, Butterfly, Moth and Bird

## 05.GLMM_visitation_indices

The script contains all the code used to run glmm’s for the three
visitation indices, which were used to assess the effect of elevation of
them.

The script contains paths to saved models in order not to have to rerun
them again.

## 06.Setup_for_index_comparisons

The script contains all the code necessary for setup of running bayesian
GLMM’s comparing the hand-pollination experiment with the visitation
data.

As mentioned above, to compare these indices, it was necessary to
aggregate the data by species and elevation level. This was because the
hand-pollination data and the visitation data were collected on
different plant individuals and so the only way to compare them was by
using bayesian statistics with measurement error and weights for each
value. This is done in the script

## 07.Comparing_indices

The script

The script contains paths to saved models in order not to have to rerun
them again.

## 08.Making_figs_main_text

The script serves to create figures found in the main text of the
manuscript

## 09.Making_tables_main_text

The script serves to create tables found in the main text of the
manuscript

## 10.Making_figs_supplementary

The script serves to create figures found in the supplementary text of
the manuscript

## 11.Making_tables_supplementary

The script serves to create tables found in the supplementary text of
the manuscript
