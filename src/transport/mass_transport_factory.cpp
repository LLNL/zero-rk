#include "stdio.h" // printf

#include "mass_transport_factory.h"
#include "constant_lewis.h"
#include "mix_avg.h"
#include "mix_avg_soret.h"

namespace transport
{

MassTransportInterface *InterfaceFactory::CreateMassBased(
    const std::string &type)
{
  if(type == "ConstantLewis") {
    return new ConstantLewis();
  } else if(type == "MixAvg") {
    return new MixAvg();
  } else if(type == "MixAvgSoret") {
    return new MixAvgSoret();
  } else {
    printf("# ERROR: In transport::InterfaceFactory::CreateMassBased(type),\n"
           "#        type = %s not supported. Returning NULL pointer.\n",
           type.c_str());
    return NULL;
  }
}

} // namespace transport
