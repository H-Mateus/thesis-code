# have a look at the data from the guide - see if the formatting is different
# and that causing problems - check the mzml file and the mzid file

library("MSnbase")
library("MSqRob")
library("msdata")
library("limma")
library("tidyverse")

# raw data
basename(fl3 <- msdata::proteomics(full.name = TRUE, pattern = "MS3TMT11"))

(rw3 <- readMSData(fl3, mode = "onDisk"))
  
# identification results
basename(idf <- msdata::ident(full.name = TRUE))

iddf <- readMzIdData(idf)
names(iddf)


# below is the head of the example mzid file
<?xml version="1.0" encoding="UTF-8"?>
  <MzIdentML id="MS-GF+" version="1.1.0" xmlns="http://psidev.info/psi/pi/mzIdentML/1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psidev.info/psi/pi/mzIdentML/1.1 http://www.psidev.info/files/mzIdentML1.1.0.xsd" creationDate="2016-10-10T13:10:37" >
  <cvList xmlns="http://psidev.info/psi/pi/mzIdentML/1.1">
  <cv id="PSI-MS" uri="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="3.30.0" fullName="PSI-MS"/>
  <cv id="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo" fullName="UNIMOD"/>
  <cv id="UO" uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="UNIT-ONTOLOGY"/>
  </cvList>
  <AnalysisSoftwareList xmlns="http://psidev.info/psi/pi/mzIdentML/1.1">
  <AnalysisSoftware version="Beta (v10072)" name="MS-GF+" id="ID_software">
  <SoftwareName>
  
# below is my mzid file  
  <?xml version="1.0" encoding="UTF-8"?>
  <MzIdentML id="MS-GF+" version="1.1.0" xmlns="http://psidev.info/psi/pi/mzIdentML/1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psidev.info/psi/pi/mzIdentML/1.1 http://www.psidev.info/files/mzIdentML1.1.0.xsd" creationDate="2020-03-26T12:59:49" >
  <cvList xmlns="http://psidev.info/psi/pi/mzIdentML/1.1">
  <cv id="PSI-MS" uri="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="3.30.0" fullName="PSI-MS"/>
  <cv id="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo" fullName="UNIMOD"/>
  <cv id="UO" uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" fullName="UNIT-ONTOLOGY"/>
  </cvList>
  <AnalysisSoftwareList xmlns="http://psidev.info/psi/pi/mzIdentML/1.1">
  <AnalysisSoftware version="Beta (v10072)" name="MS-GF+" id="ID_software">
  <SoftwareName>
  
  
# below is tail of the example
  <cvParam accession="MS:1002055" cvRef="PSI-MS" value="1.0" name="MS-GF:PepQValue"/>
  <userParam value="0" name="IsotopeError"/>
  <userParam value="HCD" name="AssumedDissociationMethod"/>
  </SpectrumIdentificationItem>
  <cvParam accession="MS:1001115" cvRef="PSI-MS" value="3017" name="scan number(s)"/>
  </SpectrumIdentificationResult>
  </SpectrumIdentificationList>
  </AnalysisData>
  </DataCollection>
  </MzIdentML>
  
# below is my file
  <cvParam accession="MS:1002055" cvRef="PSI-MS" value="1.0" name="MS-GF:PepQValue"/>
  <userParam value="0" name="IsotopeError"/>
  <userParam value="HCD" name="AssumedDissociationMethod"/>
  </SpectrumIdentificationItem>
  <cvParam accession="MS:1001115" cvRef="PSI-MS" value="3017" name="scan number(s)"/>
  </SpectrumIdentificationResult>
  </SpectrumIdentificationList>
  </AnalysisData>
  </DataCollection>
  </MzIdentML>
  
  
  # example file with line number
  105731                  <SpectrumIdentificationItem passThreshold="true" rank="1" peptide_ref="Pep4931" calculatedMassToCharge="483.8209228515625" experimentalMassToCharge="483.8291931152344" chargeState="2" id="SII_3017_1">
  105732                      <PeptideEvidenceRef peptideEvidence_ref="PepEv_1577074_4931_56"/>
  105733                      <cvParam accession="MS:1002049" cvRef="PSI-MS" value="-133" name="MS-GF:RawScore"/>
  105734                      <cvParam accession="MS:1002050" cvRef="PSI-MS" value="64" name="MS-GF:DeNovoScore"/>
  105735                      <cvParam accession="MS:1002052" cvRef="PSI-MS" value="0.010903058" name="MS-GF:SpecEValue"/>
  105736                      <cvParam accession="MS:1002053" cvRef="PSI-MS" value="31157.906" name="MS-GF:EValue"/>
  105737                      <cvParam accession="MS:1002054" cvRef="PSI-MS" value="1.0" name="MS-GF:QValue"/>
  105738                      <cvParam accession="MS:1002055" cvRef="PSI-MS" value="1.0" name="MS-GF:PepQValue"/>
  105739                      <userParam value="0" name="IsotopeError"/>
  105740                      <userParam value="HCD" name="AssumedDissociationMethod"/>
  105741                  </SpectrumIdentificationItem>
  
  # from my file
  469716                  <SpectrumIdentificationItem passThreshold="true" rank="1" peptide_ref="Pep22532" calculatedMassToCharge="575.1887817382812" experimentalMassToCharge="575.1926879882812" chargeState="2" id="SII_29706_1">
  469717                      <PeptideEvidenceRef peptideEvidence_ref="PepEv_2837391_22532_240"/>
  469718                      <cvParam accession="MS:1002049" cvRef="PSI-MS" value="-106" name="MS-GF:RawScore"/>
  469719                      <cvParam accession="MS:1002050" cvRef="PSI-MS" value="1" name="MS-GF:DeNovoScore"/>
  469720                      <cvParam accession="MS:1002052" cvRef="PSI-MS" value="0.014550735" name="MS-GF:SpecEValue"/>
  469721                      <cvParam accession="MS:1002053" cvRef="PSI-MS" value="311390.28" name="MS-GF:EValue"/>
  469722                      <cvParam accession="MS:1002054" cvRef="PSI-MS" value="1.0" name="MS-GF:QValue"/>
  469723                      <cvParam accession="MS:1002055" cvRef="PSI-MS" value="0.9787427" name="MS-GF:PepQValue"/>
  469724                      <userParam value="0" name="IsotopeError"/>
  469725                      <userParam value="HCD" name="AssumedDissociationMethod"/>
  469726                  </SpectrumIdentificationItem>
  
  
  
  ### below is a look at the mzml files from the example and my own data
  # example head (20 lines)
  <?xml version="1.0" encoding="utf-8"?>
  <indexedmzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">
  <mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" id="SP_newer_software" version="1.1.0">
  <cvList count="2">
  <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="3.79.0" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>
  <cv id="UO" fullName="Unit Ontology" version="12:10:2011" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"/>
  </cvList>
  <fileDescription>
  <fileContent>
  <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
  <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
  </fileContent>
  <sourceFileList count="2">
  <sourceFile id="RAW1" name="SP_newer_software.raw" location="file:///S:\ARON\00_PROJECTS\ARON_33_figuring_out_MS2_MS3_relation\to_send">
  <cvParam cvRef="MS" accession="MS:1000768" name="Thermo nativeID format" value=""/>
  <cvParam cvRef="MS" accession="MS:1000563" name="Thermo RAW format" value=""/>
  <cvParam cvRef="MS" accession="MS:1000569" name="SHA-1" value="1ac76961fc06ef82a79563965a9d8684f64e93d8"/>
  </sourceFile>
  <sourceFile id="SP_newer_software.mzML" name="SP_newer_software.mzML" location="file:///">
  <cvParam cvRef="MS" accession="MS:1000569" name="SHA-1" value="eac0c358fb4295b6b0a1db12d40be714c0a39a9d"/>
  
  # my file 
  <?xml version="1.0" encoding="utf-8"?>
  <indexedmzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">
  <mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" id="Mateus_iTRAQ_run1_Fr03_4ul" version="1.1.0">
  <cvList count="2">
  <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.30" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>
  <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>
  </cvList>
  <fileDescription>
  <fileContent>
  <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
  <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>
  </fileContent>
  <sourceFileList count="2">
  <sourceFile id="WIFF" name="Mateus_iTRAQ_run1_Fr03_4ul.wiff" location="file://.">
  <cvParam cvRef="MS" accession="MS:1000770" name="WIFF nativeID format" value=""/>
  <cvParam cvRef="MS" accession="MS:1000562" name="ABI WIFF format" value=""/>
  <cvParam cvRef="MS" accession="MS:1000569" name="SHA-1" value="ce00a4a8962c3b8a134bc91a0cd43af14bf59cad"/>
  </sourceFile>
  <sourceFile id="WIFFSCAN" name="Mateus_iTRAQ_run1_Fr03_4ul.wiff.scan" location="file://.">
  <cvParam cvRef="MS" accession="MS:1000770" name="WIFF nativeID format" value=""/>
  
  
  # example tail (20 lines)
  <offset idRef="controllerType=0 controllerNumber=1 scan=22927">14429957</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22928">14437254</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22929">14453495</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22930">14468105</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22931">14475082</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22932">14489375</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22933">14497603</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22934">14504576</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22935">14513014</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22936">14522351</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22937">14529858</offset>
  <offset idRef="controllerType=0 controllerNumber=1 scan=22938">14539132</offset>
  </index>
  <index name="chromatogram">
  <offset idRef="TIC">14547158</offset>
  </index>
  </indexList>
  <indexListOffset>16499467</indexListOffset>
  <fileChecksum>8bf341399c33399136f3f593d45f10da1ee66017</fileChecksum>
  </indexedmzML>
  
  # my file
  <offset idRef="sample=1 period=1 cycle=12375 experiment=1">12962634604</offset>
  <offset idRef="sample=1 period=1 cycle=12376 experiment=1">12962919675</offset>
  <offset idRef="sample=1 period=1 cycle=12377 experiment=1">12963215863</offset>
  <offset idRef="sample=1 period=1 cycle=12378 experiment=1">12963523550</offset>
  <offset idRef="sample=1 period=1 cycle=12379 experiment=1">12963813209</offset>
  </index>
  <index name="chromatogram">
  <offset idRef="TIC">12964113270</offset>
  <offset idRef="BPC">12964748476</offset>
  <offset idRef="Column Pressure (channel 1)">12965383673</offset>
  <offset idRef="Pump A Flowrate (channel 2)">12965403723</offset>
  <offset idRef="Pump B Flowrate (channel 3)">12965423790</offset>
  <offset idRef="Column Pressure (channel 4)">12965443857</offset>
  <offset idRef="Pump A Flowrate (channel 5)">12965728802</offset>
  <offset idRef="Pump B Flowrate (channel 6)">12966013764</offset>
  </index>
  </indexList>
  <indexListOffset>12966298767</indexListOffset>
  <fileChecksum>36d89c6bb3860156a3c3be32ca2789a9c8a3655a</fileChecksum>
  </indexedmzML>
  
  
  # lines 4000 - 4020 in example
  <scan>
  <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="45.56804460265" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
  <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="FTMS + p NSI Full ms [375.0000-1500.0000]"/>
  <cvParam cvRef="MS" accession="MS:1000927" name="ion injection time" value="0.043329805136" unitCvRef="UO" unitAccession="UO:0000028" unitName="millisecond"/>
  <scanWindowList count="1">
  <scanWindow>
  <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="375.0" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
  <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1500.0" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
  </scanWindow>
  </scanWindowList>
  </scan>
  </scanList>
  <binaryDataArrayList count="2">
  <binaryDataArray encodedLength="125688">
  <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
  <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
  <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
  <binary>gTyb/1o0d0DQl33uXjR3QB/zYN1iNHdAb05FzGvR45j3qsl0BnS2OqgqyXQA==</binary> #note: I deleted most of this binary thing as it was really long
</binaryDataArray>
  <binaryDataArray encodedLength="62844">
  <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
  
  # lines 4000 - 4020 in my file
  <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
  <scan>
  <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="0.32205" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
  <cvParam cvRef="MS" accession="MS:1000616" name="preset scan configuration" value="1"/>
  <scanWindowList count="1">
  <scanWindow>
  <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="400.0" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
  <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1250.0" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
  </scanWindow>
  </scanWindowList>
  </scan>
  </scanList>
  <binaryDataArrayList count="2">
  <binaryDataArray encodedLength="772192">
  <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
  <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
  <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
  # note: I skiped with line because it was way too long
  </binaryDataArray>
  <binaryDataArray encodedLength="386096">
  <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>
  
  
  