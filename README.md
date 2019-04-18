# NtProd - NtAnalyzer Workflow

README for the tHq2017 branch, describing the basic steps to run the production chain from FlatTrees to ntuples ready for MVA analysis.

*Do not forget to source :*
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh
```

*To create proxy :*
```
voms-proxy-init -voms cms -hours 192
```

## FlatTreeProducer

(( Follow instructions from [IPHCFlatTree's README](https://github.com/IPHC/IPHCFlatTree/tree/master) ))


## IPHCNtuple

Ntuple Production & Analysis codes

### Installation

```
cd /home-pbs/username
mkdir MyAnalysis/; cd MyAnalysis/

cmsrel CMSSW_9_4_3
cd CMSSW_9_4_3/src/
cmsenv

# get the code from GIT
git clone -b tHq2017  https://github.com/IPHC/IPHCNtuple.git


cd IPHCNtuple/

# NtupleProducer: produce Ntuples from FlatTrees

cd NtupleProducer/
make
cd ../../

# NtupleAnalyzer: create histograms, TTrees, ASCII, etc output from Ntuples

cd ttH/NtupleAnalyzer/
make


# to commit & push to tHq branch (from git root dir.)
git add .
git commit -m "update" //Comment your modif
git status //Check modif
git push origin tHq2017
```

### Set-up & Run NtupleProducer


```
cd /home-pbs/username/MyAnalysis/CMSSW_9_4_3/src/ttH/NtupleProducer/src
```

* **NtupleProducer.cxx** - make sure paths to JEC files for data and MC are up-to-date : 

```
jesTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("XXX", "Total")));
```






----  Produce list of FlatTree files on which to run

(NB : merging of samples, e.g. different data runs, has to be done at this step, using wildcards *)


* **input.txt** - FOR INTERACTIVE RUNNING, include directly the the Flat Tree(s) path(s) in this file, e.g. : 

```
root://sbgse1.in2p3.fr//dpm/in2p3.fr/home/cms/phedex/store/user/XXX/output_*.root
```


* **split_into_lists.zsh** - Will automatically prepare lists of FlatTree files on which to run, based on content of FlatTree dir. and maximum nof files per job :

```
...
export x509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066
...
fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/FlatTree/output_dir/"
...
```

Then execute the script to create lists of paths to FlatTree files, based on the content of your FlatTree production output directory.
By default, it will create the lists of root files in directory "lists/" : 
```
./split_into_lists.zsh
```
Then you could e.g. copy the lists for samples you're interested in in a new dir. "lists_priority", and modify run_batch.zsh accordingly.

((NB : could also list the FlatTree files yourself, and change the directory in which to look for in run_batch.zsh accordingly (cf. above) ))



```
cd /home-pbs/username/MyAnalysis/CMSSW_8_0_20/src/ttH/NtupleProducer/test
```



* **single_batch_job.zsh** - this is the code taking care of submitting a single job : 

```
...
export X509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066
...
cd /home-pbs/ntonon/tHq/CMSSW_9_4_3/src/
...
```


* **Move_File_ToGrid.sh** - after a job has been completed, this script is called to move the output to the grid (to avoid storing Ntuples on scratch1, which is an un-safe disk not meant for storing) - modify  : 

```
...
gfal-mkdir srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleProducer/$2 #Grid path where to copy outputs
gfal-copy -f $1 srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleProducer/$2/. #Grid path where to copy outputs
...
```



* **run_batch.zsh** - idem : 

```
...
cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy
...
dout="/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleProducer/test" #local dir 
dout_f="/opt/sbg/scratch1/cms/ntonon/ntuples_prod_walrus_patch2/" #output dir. (NB : only scratch1 is writable by local batch)
...
fdir=$(ls -d lists_priority*) //Can modify dir. containing list of files on which to run (e.g. to run on sub-list, ...)
...
```


### Run

*Interactively*

```
./run.zsh

//Or calling the executable yourself, e.g. : 
./NtupleProducer --file input.txt  --outfile output --tree FlatTree/tree --nmax -1 --isdata 0
```

*Launch jobs*

```
./run_batch.zsh [prod_name]
```


### Set-up & Run NtupleAnalyzer

NB : paths to weight files may need to be changed (look for ".root" in src/TTbarHiggsMultileptonAnalysis.cxx)

```
cd /home-pbs/username/MyAnalysis/CMSSW_9_4_3/src/ttH/NtupleAnalyzer/test
```
* **table.txt** - add list of samples, with cross section (pb-1) and sum of weights of events from FlatTrees (SWE, to account for efficiency from skimming, see below), e.g. : 

```
...
THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1_RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM_0000       0.7927       3495652
...
```

NB : can get SWE infos from [IPHCFlatTree twiki](https://twiki.cern.ch/twiki/bin/view/CMS/IPHCFlatTreeProduction) ) if using exact same FlatTrees.

The safer option, is to re-compute the SWE of your FlatTree files yourself (it is different from the event number, because some events have weight -1). To do this, compile and execute the code **Get_SumWeightEvents_FlatTrees.cxx** . It will ask you for the complete path to the FlatTree files of your sample, and then the total number of FlatTree files for this sample. Example : 

```
./a.out #execute the compiled code

#Enter complete path to files
> /dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/FlatTree/Walrus-patch2/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_v1_MINIAODSIM/180126_150054/0000

#Enter number of files
> 98
```



* **split_into_lists.zsh** - modify path of directory containing NtupleProducer output files, e.g. : 

```
...
fpath="/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleProducer/" #Here, the NTProducer files were stored on the grid
...
```
then run the script to produce automatically the lists of files on which to run, based on directory content and maximum nof files per job : 
```
./split_into_lists.zsh
```


* **single_batch_job.sh** - update : 

```
...
export X509_USER_PROXY=/home-pbs/ntonon/proxy/x509up_u8066

cd /home-pbs/ntonon/tHq/CMSSW_9_4_3/src/
...
```


* **Move_File_ToGrid.sh** - after a job has been completed, this script is called to move the output to the grid (to avoid storing Ntuples on scratch1, which is an un-safe disk not meant for storing) - modify  : 

```
...
gfal-mkdir srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleProducer/$2 #Grid path where to copy outputs
gfal-copy -f $1 srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleProducer/$2/. #Grid path where to copy outputs
...
```


* **run_batch.zsh** - update lumi, proxy, username, directories, e.g. : 

```
...
lumi=35900
...
cp /tmp/x509up_u8066 /home-pbs/ntonon/proxy
...
dout="/home-pbs/ntonon/tHq/CMSSW_8_0_20/src/ttH/NtupleAnalyzer/test"
dout_f="/opt/sbg/scratch1/cms/ntonon/Analyzer_ntuples_prod_walrus_patch2"
...
fdir=$(ls -d lists_NameOfYourList) //Name of the directory containing the list of files to process
...
```

-- For interactive running : 

* **run.zsh** - update lumi, xsec, nof weights, ... 

* **input.txt** - add path of input files, from NtupleProducer (for interactive running) : 


### Run

*Interactively*

```
./run.zsh
```

*Launch jobs*

```
./run_batch.zsh [prod_name]
```
