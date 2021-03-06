DIR_BAMS=/mnt/rad51/raw/human_sequence_prostate_WCDT/illumina
EXTRACT=/notebook/human_sequence_prostate_WCDT/WCDT/scripts/extract_fusion.sh
EXTRACT_DNA=/notebook/human_sequence_prostate_WCDT/WCDT/scripts/extract_DNA_translocations.sh


# 24 definitely has 1 copy loss from LOH
REGION_1=chr13:32310479-32404672
REGION_2=chr15:45000000-45000100
bash $EXTRACT_DNA -s=DTB-064-BL -d=${DIR_BAMS}/BRCA2 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-071-BL -d=${DIR_BAMS}/BRCA2 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-090-PRO -d=${DIR_BAMS}/BRCA2 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-091-BL -d=${DIR_BAMS}/BRCA2 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-165-PRO -d=${DIR_BAMS}/BRCA2 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-251-PRO -d=${DIR_BAMS}/BRCA2 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

REGION_1=chr17:43044294-43106492 
REGION_2=chr15:45000000-45000100
bash $EXTRACT_DNA -s=DTB-091-BL -d=${DIR_BAMS}/BRCA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-165-PRO -d=${DIR_BAMS}/BRCA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

# Used to generate FOXA1-intergenic-ETV1
B1=DTB-011-BL_RNA_aligned_tumor
B2=DTB-024-PRO_RNA_aligned_tumor
B3=DTB-119-PRO_RNA_aligned_tumor
REGION_1=chr14:37585500-37591000
REGION_2=chr7:13910000-13991900
bash $EXTRACT -t=${DIR_BAMS}/${B1} -d=${DIR_BAMS}/FOXA1 -o=${B1} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT -t=${DIR_BAMS}/${B2} -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT -t=${DIR_BAMS}/${B3} -d=${DIR_BAMS}/FOXA1 -o=${B3} -r1=$REGION_1 -r2=$REGION_2



# Used for SLC30
B1=DTB-055-PRO_RNA_aligned_tumor
REGION_1=chr15:45415692-45615692
REGION_2=chr7:13838878-14038878
bash $EXTRACT -t=${DIR_BAMS}/${B1} -d=${DIR_BAMS}/FOXA1 -o=${B1} -r1=$REGION_1 -r2=$REGION_2

# used for SCHLAP1 PIK3CA
B1=DTB-234-BL_RNA_aligned_tumor
REGION_1=chr2:180725105-180735105
REGION_2=chr3:179185711-179197711
bash $EXTRACT -t=${DIR_BAMS}/${B1} -d=${DIR_BAMS}/FOXA1 -o=${B1} -r1=$REGION_1 -r2=$REGION_2

# used for ACPP MYC
B1=DTB-005-BL_RNA_aligned_tumor
B2=DTB-035-BL_RNA_aligned_tumor
REGION_1=chr3:132321189-132359962
REGION_2=chr8:127724991-127744991
bash $EXTRACT -t=${DIR_BAMS}/${B1} -d=${DIR_BAMS}/FOXA1 -o=${B1} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT -t=${DIR_BAMS}/${B2} -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-035-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

# VRK3 MYC
B1=DTB-205-BL_RNA_aligned_tumor
REGION_1=chr19:50003926-50023926
REGION_2=chr8:127633147-127733147
bash $EXTRACT -t=${DIR_BAMS}/${B1} -d=${DIR_BAMS}/FOXA1 -o=${B1} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-205-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

#DTB-151-BL	SCHLAP1	FOXA1	chr2	180895623	chr14	37574528
B1=DTB-151-BL_RNA_aligned_tumor
REGION_1=chr2:180885623-180905623
REGION_2=chr14:37474528-37674528 
bash $EXTRACT -t=${DIR_BAMS}/${B1} -d=${DIR_BAMS}/FOXA1 -o=${B1} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-151-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

REGION_1=chr21:41490359-41510359
REGION_2=chr21:38482410-38502410  
bash $EXTRACT_DNA -s=DTB-083-BL -d=${DIR_BAMS}/FOXA1 -o=DTB-083-BL_TE_FUSION -r1=$REGION_1 -r2=$REGION_2

REGION_1=chr8:127137709-127157709
REGION_2=chr15:45509237-45529237   
bash $EXTRACT_DNA -s=DTB-042-BL -d=${DIR_BAMS}/FOXA1 -o=DTB-083-BL_TE_FUSION -r1=$REGION_1 -r2=$REGION_2

# 24 definitely has 1 copy loss from LOH
REGION_1=chr13:32315479-32399672
REGION_2=chr15:45000000-45000100
bash $EXTRACT_DNA -s=DTB-024-PRO -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2


8q24 region


REGION_1=chr8:127400510-127473443
REGION_2=chr8:127012153-127022014
bash $EXTRACT_DNA -s=DTB-018-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-022-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-023-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-040-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-055-PRO -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-060-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-074-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-080-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

bash $EXTRACT_DNA -s=DTB-067-PRO -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-098-PRO2 -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-104-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-126-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-132-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-143-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-146-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-165-PRO -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-170-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-173-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-188-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-194-PRO -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-206-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

bash $EXTRACT_DNA -s=DTB-188-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-213-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-251-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2



bash $EXTRACT_DNA -s=DTB-063-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2
bash $EXTRACT_DNA -s=DTB-060-BL -d=${DIR_BAMS}/FOXA1 -o=${B2} -r1=$REGION_1 -r2=$REGION_2


LP6008330-DNA-A01_WFR2351179_S1.genome.vcf.gz 

[1] "DTB-035-BL"  "DTB-059-BL"  "DTB-063-BL"  "DTB-100-BL"  "DTB-101-BL"  "DTB-112-BL"  "DTB-119-PRO"
 [8] "DTB-127-PRO" "DTB-129-BL"  "DTB-170-BL"  "DTB-173-BL"  "DTB-176-BL"  "DTB-183-BL"  "DTB-190-BL" 
[15] "DTB-193-BL"  "DTB-205-BL"  "DTB-214-BL"  "DTB-216-PRO" "DTB-220-BL"  "DTB-232-PRO" "DTB-252-BL" 
[22] "DTB-260-BL"  "DTB-265-PRO"




# CDK12

 DTB-214-BL
 DTB-183-BL

DTB-063-BL "chr17:39494625_C>T_stop_gained_35_60" 

REGION_1=chr17:38000000-40000000
REGION_2=chr7:55019031-55020031
bash $EXTRACT_DNA -s=DTB-214-BL -d=${DIR_BAMS}/CDK12 -o=${B2} -r1=$REGION_1 -r2=$REGION_2

REGION_1=chr17:38000000-40000000
REGION_2=chr7:55019031-55020031
bash $EXTRACT_DNA -s=DTB-183-BL -d=${DIR_BAMS}/CDK12 -o=${B2} -r1=$REGION_1 -r2=$REGION_2


# checked DTB-214-BL sequence manually. Only identified the single mutation, no evidence for LOH or copy number alteration.
# SNPs that were near 50% VAF in germline also near 50% in tumor.

DTB-214-BL	CDK12	missense_variant	chr17	39492849	G	A	93	37
DTB-214-BL	CDK12	missense_variant	chr17	39511622	G	C	102	26
# DTB-063-BL: stop gain CDK12 and copy loss
# DTB-214-BL: two mutations in CDK12
# DTB-183-BL: LOH and single mutations

DTB-003-BL chr21 38498001 41506000 LOSS 0.6958763      +
REGION_1=chr21:38488001-38508001
REGION_2=chr21:41496000-41516000   
bash $EXTRACT_DNA -s=DTB-003-BL -d=${DIR_BAMS}/FOXA1 -o=DTB-003-BL_TE_FUSION -r1=$REGION_1 -r2=$REGION_2

Read name = 943:330:HCJV3CCXY:3:2218:11820:72948
Sample = LP6008328-DNA-A01_WFR2307750
Read group = 2
----------------------
Location = chr21:38,498,139
Alignment start = 38,498,116 (-)
Cigar = 42M108S
Mapped = yes
Mapping quality = 60
Secondary = no
Supplementary = no
Duplicate = no
Failed QC = no
----------------------
Base = T
Base phred quality = 41
----------------------
Mate is mapped = yes
Mate start = chr21:38497800 (+)
Insert size = -149
Second in pair
Pair orientation = F1R2
----------------------
SA = chr12,71004887,-,42S108M,60,0;
BC = none
OC = 42M11C32506729D108M
RG = 2
NM = 0
SM = 65535
AS = 1055
-------------------
Alignment start position = chr21:38498116
AATCCACAGCACTCAACAGTGAATTTAGAAACCCCAATTTTTGATGTAATCAAGTTAAAATGAGGTCATAATGGATTCAGTCGGGACCTAATCCAATGATCCATGCCCTAATAAAAAGAGAGAAATTTGTGCTTTTTAAAAACAAATGTT

Read name = 943:330:HCJV3CCXY:3:2112:21663:69572
Sample = LP6008328-DNA-A01_WFR2307750
Read group = 2-6170EA33
----------------------
Location = chr21:41,506,796
Alignment start = 41,506,796 (+)
Cigar = 132S18M
Mapped = yes
Mapping quality = 60
Secondary = no
Supplementary = yes
Duplicate = no
Failed QC = no
----------------------
Base = C
Base phred quality = 37
----------------------
Mate is mapped = yes
Mate start = chr21:41506865 (-)
Insert size = -29501819
First in pair
Pair orientation = R2F1
----------------------
SA = chr12,71008835,+,132M18S,60,0;
BC = none
OC = 132M20C29502171B18M
RG = 2-6170EA33
NM = 0
SM = 65535
AS = 358
-------------------
Alignment start position = chr21:41506796
TGCAGTTCAGTTTACACTGCAACTGAAACTTGGCTTCTAGCTTCTTGTTGACAATGTGCTCATATAGGCAAAAAGCTTTCCTTAATACATTTTGTAATACAATTCAGAAACATTAAAAAATAAGTGTTTATC CTGCGTAGCCTTTGTAAT