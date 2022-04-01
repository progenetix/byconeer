BASEDIR=$(dirname $0)
UNITYPATH=$BASEDIR/../bycon/beacon_schemas

# initial conversion from separate schemas
BEACONMODELPATH=/Users/$USER/Github/beacon-v2/models
BEACONFRAMEWORKPATH=/Users/$USER/Github/beacon-v2/framework

for KIND in src json
do
	mkdir -p $UNITYPATH/models/$KIND
	mkdir -p $UNITYPATH/framework/$KIND	
done
