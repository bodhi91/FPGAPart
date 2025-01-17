#pragma once

namespace par {
class PartitionMgr;
}

namespace ord {

class OpenRoad;

par::PartitionMgr* makePartitionMgr();

void initPartitionMgr(OpenRoad* openroad);

void deletePartitionMgr(par::PartitionMgr* partitionmgr);

}  // namespace ord
