#Need to keep this file in ./test, needed to move outputs of batch jobs to
#correct location

if [ "$1" == "" ]
then
echo "Missing arg : file name"
elif [ "$2" == "" ]
then
echo "Missing arg : destination"
elif [ "$3" == "" ]
then
echo "Missing arg : version"
fi

#Copy file from scratch1 to grid (overwrite)
gfal-mkdir srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/$3
gfal-mkdir srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/$3/$2
gfal-copy -f $1 srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/$3/$2/.

#Then remove file from scratch1
rm $1
