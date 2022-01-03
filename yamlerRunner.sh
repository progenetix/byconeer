BYCONPATH=~/groupbox/dbtools/bycon
BYCONEERPATH=~/groupbox/dbtools/bycon
BYCONRSRCPATH=~/groupbox/dbtools/bycon/beaconServer/rsrc
BEACONMODELPATH=~/GitHub/beacon-v2-Models/BEACON-V2-draft4-Model
BEACONFRAMEWORKPATH=~/GitHub/beacon-framework-v2

./yamler.py -i $BEACONMODELPATH -t json -o $BYCONRSRCPATH/beacon-v2-dev/model
./yamler.py -i $BEACONFRAMEWORKPATH -t json -o $BYCONRSRCPATH/beacon-v2-dev/framework
