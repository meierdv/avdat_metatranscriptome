/*looking for RuBisCO*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Ribulose bisphosphate carboxylase%'

/*looking for ATP-citrate lyase (rTCA cycle key enzyme)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%2.3.3.8%' OR
RAST like '%ATP-citrate%lyase%' OR
Pfam like '%ATP citrate lyase citrate-binding%'

/*searching for 2-oxoglutarate oxidoreductase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%1.2.7.3%' OR
Pfam like '%Pyruvate:ferredoxin oxidoreductase core domain%'

/*searching for fumarate reductase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%1.3.5.4%' OR
RAST like '%fumarate%reductase%'

/*looking for CO-dehydrogenase/acetylCoA-decarboxylase fusion enzyme, reductive Acetyl-coA carbon fixation key enzyme*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%1.2.7.4%' OR
RAST like '%acetyl%decarboxylase%' OR
RAST like '%carbon%dehydrogenase%'

/*searching biotin carboxylase* (3-hydroxypropanoate cycle, CO2-fixation in Chloroflexi, many multi-functional enzymes all needed, no indicator key enzyme)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%6.3.4.14%'

/*sarching malonyl-CoA reductase (3-hydroxypropanoate cycle, CO2-fixation in Chloroflexi)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%1.2.1.75%' OR
RAST like '%1.1.1.298%' OR
RAST like '%molon%reductase%' OR
RAST like '%3-hydroxypropionate%oxidoreductase%' OR
RAST like '%3-hydroxypropionate%dehydrogenase%'

/*searching propionyl-CoA synthase (3-hydroxypropanoate cycle, CO2-fixation in Chloroflexi)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%6.2.1.36%' OR
RAST like '%4.2.1.116%' OR
RAST like '%hydroxypropionyl%'

/*searching propionyl-CoA carboxylase (3-hydroxypropanoate cycle, CO2-fixation in Chloroflexi)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%6.4.1.3%' OR
RAST like '%Propionyl%carboxylase%'

/*looking for photosystem components*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%Photosystem I %'

/*looking for photosystem components*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%Photosystem II %'

/*looking for chlorophyll synthesis*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%chlorophyll%' and
RAST not like '%Photosystem II%'

/*looking for blue light sensors*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%blue%' OR
Uniprot like '%blue light%'

/*looking for carotene synthesis*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Beta-carotene%'
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%carotenoid biosynthesis%'
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%caroten%ketolas%'

/*looking for rhodopsins*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Bacteriorhodopsin%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Heliorhodopsin%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%sensory rhodopsin%'

/*looking for hydrogenases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE
RAST like '% hydrogenase%subunit%' AND Pfam like '% hydrogenase%' OR
Pfam like '%Fe% hydrogenase%' OR
Pfam like 'Nickel% hydrogenase%' OR
Pfam like '%| Nickel% hydrogenase%' OR
Uniprot like 'Hydrogenase%subunit%' AND Pfam like '% hydrogenase%' OR
Uniprot like '%Fe%hydrogenase%subunit%' AND Pfam like '% hydrogenase%'

/*looking for CO-dehydrogenase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%carbon%dehydrogenase%chain%'

/*looking for soluble Methane monooxygenase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%methane%component A%alpha%' OR
Pfam like '%Toluene hydroxylase%' AND RAST not like '%methane%component A%beta%' OR
EggNOG is '0XRJ2'


/*looking for putative (!) particulate Methane monooxygenase/Ammonia monooxygenase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%pmoa%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%ammonia monooxygenase%'


/*looking for sulfur oxidation*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%Flavocytochrome%sulfide%dehydrogenase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%sulfur%oxidation%' OR
Pfam like '%sulfur%oxidation%'

/*looking for terminal oxidases (Cytochrome C)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%Cytochrome c oxidase%1.9.3.1%'

/*looking for terminal oxidases (Cytochrome C)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%Menaquinone-cytochrome c%'

/*looking for terminal oxidases (Cytochrome bd)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Cytochrome bd terminal oxidase%'

/*looking for alternative oxidases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%alternative%oxidase%'

/*looking for nitrite reducatases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%nitrite%reductase%NAD%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%ferredoxin%nitrite%reductase%'

/*looking for lactic acid fermentation*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%lactate dehydrogenase%' OR
Pfam like '%lactate/malate dehydrogenase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%pyruvate%ferredoxin%oxidoreductase%'

/*looking for mixed acid fermentation*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%formate lyase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%formate%hydrogen%lyase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%acetate kinase%' OR
EggNOG is 'COG0282'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%phosphate acetyltransferase%' OR
EggNOG is 'COG0280'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Citrate synthase%' AND
(RAST like '%(EC 2.3.3.1%' OR
EggNOG is 'COG0372')


SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%(EC 4.2.1.3)%' OR
EggNOG is 'COG1048'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%(EC 1.1.1.42)%' OR
EggNOG is 'COG2838' OR
EggNOG is 'COG2838'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%acetaldehyde dehydrogenase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Iron-containing alcohol dehydrogenase%' AND
(RAST like '%Alcohol dehydrogenase%' OR
EggNOG is 'COG1012' OR
EggNOG is 'COG1454﻿')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Aldehyde dehydrogenase family%' AND
(RAST like '%(EC 1.2.1.3)%' OR
EggNOG is 'COG1012')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
Pfam like '%Fumarase%' AND
(RAST like '%(EC 4.1.1.31)%' OR
EggNOG is 'COG2352﻿﻿')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Phosphoenolpyruvate carboxylase%' AND
(RAST like '%(EC 4.1.1.31)%' OR
EggNOG is 'COG2352﻿')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%malate dehydrogenase%' AND
(RAST like '%(EC 1.1.1.37)%' OR
EggNOG is 'COG0039﻿')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
Pfam like '%Fumarase%' AND
(RAST like '%(EC 4.2.1.2)%' OR
EggNOG is 'COG0114﻿')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
Pfam like '%Fumarate reductase%' AND 
(RAST like '%(EC 1.3.5.4)%' OR
EggNOG is 'COG1053﻿')

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
(Pfam like '%2Fe-2S iron-sulfur cluster binding domain%' OR
Pfam like '%4Fe-4S dicluster domain%') AND
(RAST like '%(EC 1.3.5.4)%' OR
EggNOG is 'COG0479﻿')

/*looking for nitrogenases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Uniprot like '%nitrogenase%' OR
RAST like '%nitrogenase%'

/*looking for NifK homologues*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%chlorophyllide%subunit B %'

/*looking for NifD homologues*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%chlorophyllide%subunit N %'

/*looking for SusC/D-like transporters*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%TonB dependent receptor %' OR
Pfam like '%SusD%' OR
RAST like '%SusC%' OR
RAST like '%SusD%' OR
EggNOG is '0XNNV' OR
EggNOG is '0XP53' OR
EggNOG is '0XSVM' OR
EggNOG is '0YT3X'


/*looking for Carbohydrate Binding Motifs*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
HMM_vs_CAZY like '%CBM%'

/*looking for Glycosyl Hydrolases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
HMM_vs_CAZY like '%GH%'

/*looking for Proteases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
MEROPS_diamond not like ''

/*looking for sugar transporters*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%port%' AND RAST like '%ose%'

/*looking for peptide transporters*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%port%' AND RAST like '%peptide%'

/*looking for amino acid transporters*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%port%' AND RAST like '%amino%acid%'

/*looking for lactate permease*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%permease%' AND RAST like '%lactate%'
OR Pfam like '%L-lactate permease%'

/*looking for other permeases*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%permease%' AND RAST not like '%port%' AND RAST not like '%PTS%'

/*looking for PTS systems*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%PTS%' AND not RAST like '%nitroge%regulatory%'

/*looking for pilin*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%pil%' and RAST not like '%twitching%'

/*looking for twitching motility*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%twitching%'

/*looking for gliding motility*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%Gld%' OR
RAST like '%gliding%' OR
Uniprot like '%Gld%' OR
Uniprot like '%Gliding%' OR
EggNOG_function like '%Gld%' OR
EggNOG_function like '%gliding%'

/*looking for different superoxide dismutatses*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE
RAST like '%Nickel%superoxide%' OR
EggNOG like '%11IW1%' OR
Pfam like '%Nickel-containing superoxide dismutase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
EggNOG like '%COG0605%' AND
RAST like '%superoxide%[Fe]%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
EggNOG like '%COG0605%' AND
RAST like '%superoxide%[Mn]%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
EggNOG like '%COG0605%' AND
RAST like '%superoxide%[Mn/Fe]%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE
EggNOG like '%COG2032%' AND
RAST like '%superoxide%[Cu-Zn]%' OR
Pfam like '%Copper/zinc superoxide dismutase (SODC)%'

/*looking for catalase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%catalase%' OR
Pfam like '%catalase%'

/*looking for Mn-transporters*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%transporter%' AND Rast like '%Mn%Zn%' OR
RAST like '%transporter%' AND Rast like '%Manganes%'

/*looking for Cytochrome C peroxidase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%cytochrome% peroxidase%' OR
Pfam like '%cytochrome% peroxidase%'

/*looking for animal haem peroxidase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%animal% peroxidase%'

/*looking for glyoxilate cycle enzymes*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Isocitrate lyase family%' AND
RAST like 'Isocitrate lyase%'

SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Malate synthase%' AND
RAST like 'Malate synthase%'

/*Rubrerythrins*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%Rubrerythrin%'


/*looking for DNA photolyase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%photo%lyase%' OR
Pfam like '%photolyase%'

/*looking for other DNA damage repair*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%DNA%' and RAST like '%repair%'

/*looking for other DNA damage repair*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%UV%' AND
RAST like '%repair%' OR
Pfam like '%UV%endonuclease%'

/*looking for DNA polymerase PolA*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG like '%COG0749%' OR
EggNOG like '%COG0258%'

/*looking for DNA polymerase IV*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG like 'COG0389'

/*looking for Exodeoxyribonuclease VII*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG like 'COG1570'

/*searching for capsular polysaccharide synthesis*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Uniprot like '%exopolysaccharide%' OR
RAST like '%exopolysaccharide%' OR
EggNOG_function like '%exopolysaccharide%' OR
Pfam like '%Exopolysaccharide%'

/*searching for glycogen metabolism*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Uniprot like '%glycogeb%' OR
RAST like '%glycogen%' OR
EggNOG_function like '%glycogen%' OR
Pfam like '%glycogen%'

/*searching for PHAs*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Uniprot like '%glycogeb%' OR
RAST like '%PHB%' OR
RAST like '%Poly%hydroxyalkanoate%' OR
RAST like '%Poly%hydroxybutyrate%' OR
RAST like '%Polyhydroxyalkanoic%'
Pfam like '%PHB%' OR
Pfam like '%hydroxyalkanoate%' OR
Pfam like '%hydroxybutyrate%' OR
EggNOG_function like '%glycogen%' OR
Pfam like '%glycogen%'

/*searching for vitrification related compounds synthesis*/

/*searching for trehalose synthase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE
RAST like '%trehalose synthase%' OR
EggNOG is 'COG3280' OR
EggNOG is 'COG3280'

/*searching for mannosylglycerate synthase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE
CAZY_HMM like 'GT81%'

/*searching for glycine/betaine transporter*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE
RAST like '%betaine%transport%' OR
EggNOG is 'COG1292' OR
EggNOG is 'COG2113' OR
EggNOG is 'COG1174' OR
EggNOG is 'COG4176' OR
EggNOG is 'COG4175'



/*searching anything sporulation related*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Uniprot like '% spore %' OR
Uniprot like '% sporu%' OR
Pfam like '% spore %' OR
Pfam like '% sporu%' OR
RAST like '% spore %' OR
RAST like '% sporu%' OR
EggNOG_function like '% spore %' OR
EggNOG_function like '% sporu%' OR
Uniprot like 'spore %' OR
Uniprot like 'sporu%' OR
Pfam like 'spore %' OR
Pfam like 'sporu%' OR
RAST like 'spore %' OR
RAST like 'sporu%' OR
EggNOG_function like 'spore %' OR
EggNOG_function like 'sporu%'

/*searching for histone-like DNA conservation proteins*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG_function like 'integration host factor' OR
EggNOG_function like '%Lsr2-like protein%' OR
EggNOG like '%0XPFG%'

/*Assymetric cell division FstX*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG2177'

/*Streptomyces sporulation and cell division protein, SsgA*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%SsgA%'

/*DNA-binding transcriptional regulator WhiA, involved in cell division*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG1481'

/*Transcriptional regulator WhiB*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '11TPB' OR
EggNOG is '11ZPZ'

/*Transcriptional regulator WhiD*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '11U54'

/*RNA polymerase sporulation specific sigma factor SigH*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
RAST like '%RNA polymerase sporulation specific sigma factor SigH%' OR
Uniprot like '%RNA polymerase sporulation sigma factor SigH%' OR
Uniprot like '%RNA polymerase, sigma 30 subunit, SigH%'


/*Cell wall hydrolase CwlJ, involved in spore germination*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG3773'


/*Immunoglobulin-like domain of bacterial spore germination*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '0ZYB3' OR
EggNOG is '0ZVA2' OR
EggNOG is '0ZWGT' OR
EggNOG is '0ZVUJ' OR
EggNOG is '126R3'

/*Lipoprotein LpqB, GerMN domain*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '0ZF99' OR
EggNOG is '11S85' OR
EggNOG is '0Y6FI' OR
EggNOG is '11K5J'

/*Polysaccharide deacetylase family sporulation protein PdaB*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '112BD'

/*Protein of unknown function identified by role in sporulation (SpoVG)*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
Pfam like '%SpoVG%' OR
EggNOG is '11SKX'

/*Sporulation-control protein spo0M*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG4326' OR
Pfam like '%SpoOM protein%'

/*spore coat polysaccharide biosynthesis protein*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '11SMP'

/*Spore coat sialic acid synthase SpsE*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG2089'

/*Spore coat sialic acid synthase SpsF*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG1861'

/*Spore coat protein CotF*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG5577'

/*Spore coat protein CotH*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '0YUJE' OR
EggNOG is 'COG5337' OR
EggNOG is '0XP9X'

/*Spore germination protein GerM/Lipoprotein LpqB, GerMN domain*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG5401'

/*Spore germination protein YaaH*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG3858'

/*Spore germination protein CgeB*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG4641'

/*Spore maturation protein SpmB*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG0700'

/*Uncharacterized conserved protein YeaH/YhbH, required for sporulation*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG2718'

/*Sporulation-specific N-acetylmuramoyl-L-alanine amidase*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is '0XTIC'

/*Stage II sporulation protein D*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG2385'

/*Stage II sporulation protein M*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG1300'

/*Stage III sporulation protein AA*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG3854'

/*Stage V sporulation protein R*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG2719'

/*Stage V sporulation protein S*/
SELECT * from Avdat_per_bin_normalized_TPMs WHERE 
EggNOG is 'COG2359'
