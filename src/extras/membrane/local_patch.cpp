#include "pteros/extras/membrane/local_patch.h"
#include <ranges>

using namespace  std;

void pteros::LocalPatch::sort_and_remove_duplicates(){
    ranges::sort(neib_id);
    // Remove duplicates
    auto const [del1,del2] = ranges::unique(neib_id);
    neib_id.erase(del1,del2);
}
