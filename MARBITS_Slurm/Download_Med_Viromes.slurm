#!/bin/sh
#
# SLURM HEADER
#############################################
# JOB INFO
#SBATCH --account=emm4
#SBATCH --time=04-00:00:00
#SBATCH --job-name=Download_BOV_Dataset
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.out
#SBATCH --error=/mnt/lustre/scratch/fcoutinho/Job_Logs/jobLog_%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=felipehcoutinho@gmail.com

##############################################

wget https://de.cyverse.org/dl/d/28C309E5-EC1A-4953-9FE6-D4DF2BD86DE4/MSP_DeepMed_215_10312.5.156157.CGTACTA-CTCTCTA.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/C58436AA-5136-4B7F-A9B5-23E5C1AE272B/MSP_DeepMed_215_10312.5.156157.CGTACTA-CTCTCTA.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/C181A566-9F6A-4198-B784-2C3536445287/MSP_DeepMed_231_10312.6.156164.AGGCAGA-TATCCTC.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/F3B02FD5-A0CC-4280-A230-0F6A21C91F85/MSP_DeepMed_231_10312.6.156164.AGGCAGA-TATCCTC.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/E0F991B3-BD85-4641-8777-00E31925CE3F/MSP_DeepMed_51_10312.6.156164.AGGCAGA-AGAGTAG.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/107BB694-F07D-4BD2-937B-7C9B279D7D42/MSP_DeepMed_51_10312.6.156164.AGGCAGA-AGAGTAG.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/C1EE0F6B-A920-428A-9478-5508F5E622C2/MSP_DeepMed_904_10312.6.156164.CGTACTA-TATCCTC.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/6F7102B9-27AD-45D0-86D5-1485853F7478/MSP_DeepMed_904_10312.6.156164.CGTACTA-TATCCTC.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/B90B4D07-859A-4590-8F28-D425CBF24085/MSP_DeepMed_906_10359.2.158623.TAGGCAT-AGAGTAG.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/592AA58F-3A11-46A4-9B21-69E34C5FFFFF/MSP_DeepMed_906_10359.2.158623.TAGGCAT-AGAGTAG.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/6E1B6AE4-1208-4291-8F8A-056A565C5B79/MSP_DeepMed_908_10312.5.156157.TAAGGCG-CTCTCTA.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/51948231-B7B3-46F5-99E8-48C915940EFF/MSP_DeepMed_908_10312.5.156157.TAAGGCG-CTCTCTA.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/426BB058-6319-4A8F-B7C1-62CDFDA6382D/MSP_DeepMed_Geostar_10359.2.158623.GGACTCC-TATCCTC.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/62BFD995-EA2B-4232-A5F6-7740B1958D19/MSP_DeepMed_Geostar_10359.2.158623.GGACTCC-TATCCTC.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/064B9B7B-26D7-479E-A389-4A6666901EEA/MSP_DeepMed_M9_10312.6.156164.AGGCAGA-GCGTAAG.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/B270B970-6066-47C7-846B-68BDAB62A04C/MSP_DeepMed_M9_10312.6.156164.AGGCAGA-GCGTAAG.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/3A8350C2-412D-4FB8-B991-C2ABA3C2325D/MSP_DeepMed_b05_10312.5.156157.CGTACTA-AGAGTAG.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/71B87B1E-DABF-4A5E-A5B5-F9FEE0EA1CB5/MSP_DeepMed_b05_10312.5.156157.CGTACTA-AGAGTAG.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/ACCE7E9B-24C5-47E7-A10E-CBEDBE429FE8/MSP_DeepMed_s6_10312.5.156157.TAAGGCG-TATCCTC.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/A31D7A6D-2C25-4DB1-95B6-2DE22E575753/MSP_DeepMed_s6_10312.5.156157.TAAGGCG-TATCCTC.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/A0F7B70C-43DD-40BE-AC3E-92BCB967C7DE/MSP_DeepMed_s12_10359.2.158623.GGACTCC-GCGTAAG.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/EFE7BA40-E028-47FA-BCD3-9281A4EABA87/MSP_DeepMed_s12_10359.2.158623.GGACTCC-GCGTAAG.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/73325443-71C4-46DF-8519-A8A7698AFF29/MSP_DeepMed_s16_10312.6.156164.TCCTGAG-AGAGTAG.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/02F702E5-A1EF-41B6-A8C6-5E49A88F5E6E/MSP_DeepMed_s16_10312.6.156164.TCCTGAG-AGAGTAG.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/E3969EBD-B831-4B74-A61B-0B1BBE269707/MSP_DeepMed_s17_10359.2.158623.GGACTCC-TATCCTC.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/4892B6CC-6AA5-40E0-816C-63FC4D330F68/MSP_DeepMed_s17_10359.2.158623.GGACTCC-TATCCTC.anqdpht_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/24639AFF-1EF8-41D1-A0E5-6D9C10D7DA95/MSP_DeepMed_s2_10359.2.158623.TAGGCAT-CTCTCTA.anqdpht_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/69D521AF-5458-4847-A4DB-8871C9EA06AD/MSP_DeepMed_s2_10359.2.158623.TAGGCAT-CTCTCTA.anqdpht_R2_split_t_paired.fastq.gz

wget https://de.cyverse.org/dl/d/A4237DE7-176B-426F-AC43-9455D4176EC0/MSP103_7030.1.59052.GGACC_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/701D2E99-8AF5-439B-A89D-0BF8E096A2F0/MSP103_7030.1.59052.GGACC_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/73123D42-D66B-432C-B2D6-4888C760EDBA/MSP109_7036.1.59404.TTCAG_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/9E2B8F91-4652-46B1-93F2-7C2FC327259D/MSP109_7036.1.59404.TTCAG_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/131ADBD5-A63E-4FD8-87BB-CB0E2D3648CE/MSP112_7153.1.61927.TGTGA_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/E165A8D6-E0E1-4D04-BF8D-2CEC2B1A3F2C/MSP112_7153.1.61927.TGTGA_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/A0A9B2E7-182A-4FD2-BF92-8A3B21570EC5/MSP114_7151.1.61919.ACAAA_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/29EF89B1-2C4B-44D6-9086-0ABEF7C5296C/MSP114_7151.1.61919.ACAAA_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/E6DB1B58-C9DB-457B-B44F-C58C2E725F27/MSP118_7154.1.61923.CACCT_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/CEA47D41-ABDB-4B5B-8E3A-070F4D15AA4C/MSP118_7154.1.61923.CACCT_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/44BA2DFA-A1DC-49F8-970C-02333775FF93/MSP121_7152.1.61920.GTGGC_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/431CBF2E-1219-45DF-A316-07BEAC66964A/MSP121_7152.1.61920.GTGGC_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/693FF29B-3459-4054-A3ED-E7F5F61C96BA/MSP131_7037.1.59410.ACACG_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/694B74B9-5963-4E18-B888-9E9B23A00970/MSP131_7037.1.59410.ACACG_R1_split_t_unpaired.fastq.gz
wget https://de.cyverse.org/dl/d/B2335DF0-5B8E-4D01-BE48-5A5CB94479A1/MSP131_7037.1.59410.ACACG_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/79A6284C-F384-4AC1-A8D0-42F1CA1F2FE7/MSP134_7035.1.59396.CACAC_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/FEC61351-A1AE-42C6-8A5C-02C8D9DE32B9/MSP134_7035.1.59396.CACAC_R1_split_t_unpaired.fastq.gz
wget https://de.cyverse.org/dl/d/0E69775E-E64A-4467-BE7F-383EC5597CE6/MSP134_7035.1.59396.CACAC_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/1EBFCF31-39E1-400D-AEEC-FDBA950AECCC/MSP144_7034.1.59397.GTGTT_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/29E1DC26-E8FB-40A4-A3E7-F23556CD563D/MSP144_7034.1.59397.GTGTT_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/8AD34F84-70E9-4590-943E-5F2AB9A5BB42/MSP146_7058.1.60048.TGTGA_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/D6268CBB-96F6-4771-A732-F110705D9B07/MSP146_7058.1.60048.TGTGA_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/EFB9F98B-56D1-452C-AE8A-39EC4FB10FA9/MSP81_6982.1.58137.GGACC_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/70BDF269-AC10-4128-A3B4-B0DF747A87DD/MSP81_6982.1.58137.GGACC_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/D6F2FD25-1492-49FE-BEB1-3A787E5C52F9/MSP82_6985.1.58159.TTCAG_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/0A71EB03-B61A-4DDB-A2AD-139FCDCFA283/MSP82_6985.1.58159.TTCAG_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/B9F17EED-EF25-4F3F-BD70-7497448464E8/MSP91_7029.1.59050.AAGGG_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/4B994760-665C-4AF7-B7E9-7B76DB07E78E/MSP91_7029.1.59050.AAGGG_R2_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/35135EB8-6056-442C-8F3D-1FAF814E0875/MSP97_7028.1.59051.CCTTC_R1_split_t_paired.fastq.gz
wget https://de.cyverse.org/dl/d/F4D8A8FE-E811-49DE-ACF6-9518A50AD971/MSP97_7028.1.59051.CCTTC_R2_split_t_paired.fastq.gz