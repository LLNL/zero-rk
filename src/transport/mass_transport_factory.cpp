#include "stdio.h" // printf

#include "mass_transport_factory.h"
#include "constant_lewis.h"
#include "mix_avg.h"
#include "mix_avg_soret.h"
#include "flexible_transport.h"

namespace transport
{

MassTransportInterface *InterfaceFactory::CreateMassBased(
    const std::string &type)
{
  if(type == "ConstantLewisOld") {
    return new ConstantLewis();
  } else if(type == "MixAvgOld") {
    return new MixAvg();
  } else if(type == "MixAvgSoretOld") {
    return new MixAvgSoret();
  } else if(type == "ConstantLewis") {
    FlexibleTransport* ptr = new FlexibleTransport();
    return ptr;
  } else if(type == "MixAvg") {
    FlexibleTransport* ptr = new FlexibleTransport();
    ptr->SetMixAvg(true);
    return ptr;
  } else if(type == "MixAvgSoret") {
    FlexibleTransport* ptr = new FlexibleTransport();
    ptr->SetMixAvg(true);
    ptr->SetSoret(true);
    return ptr;
  } else if(type == "Flexible") {
    return new FlexibleTransport();
  } else {
    printf("# ERROR: In transport::InterfaceFactory::CreateMassBased(type),\n"
           "#        type = %s not supported. Returning NULL pointer.\n",
           type.c_str());
    return NULL;
  }
}

} // namespace transport
