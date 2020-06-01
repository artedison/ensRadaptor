# explanation of provided data
## the latest ens code
./enscode

## collected kinetic information for Neurospora crassa
./kin_infor

If the user need to prepare this for other specific organisms, they need parepare a list of km and kcat range and specific value in databases just as all.kine.RData or Neurospora crassa.kine.new.RData. The sublist will be with name of reaction.

## the stored measurement files
./measurements

If the user need to prepare this for input, they need to parepare quantification of each compounds as a dataframe in a list. Multiple experiment or replicate can be included as in "exp".

## the pathway reaction list
./pathway

The format should be easy to follow and names of reaction, compound, and enzymes should in agreement with the script. The start of each line "Reaction (\d)" isn't read.

## template format under ./template_format
This is used by function to format output

spec.metabolite.tab: the default metabolite formatting template under ens.i01

spec.enz.massscal.tab: the enzyme formatting template under ens.i01 with supporting on biomass proportions (all enzyme concentration are proportional with biomass)

reac.rev.mm.tab: the reversible reaction formatting template under ens.i01

reac.irrev.mm.tab: the irreversible reaction formatting template under ens.i01

comb_glycolysis.tab: the combined/simplified reaction for glycolysis ens.i01

addon.comb.meas.tab: the template for formulating combined measurement: the measured signal is a weighted sum of internal and external concentrations

addon.scalenz.tab: the template for formualting biomass proportions for all enzymes

addon.aero.vs.anaero.tab: the template for formualting different reactions under aerobic and anaerobic conditions
