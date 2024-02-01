#include <ost/log.hh>
#include <ost/conop/conop.hh>
#include <ost/io/mol/entity_io_mae_handler.hh>

using namespace ost;

int main(int argc, char** argv)
{
  if(argc>1) {
    mol::EntityHandle eh=io::LoadMAE(argv[1]);
  }

  return 0;
}
