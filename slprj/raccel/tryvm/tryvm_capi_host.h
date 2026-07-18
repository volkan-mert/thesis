#ifndef tryvm_cap_host_h__
#define tryvm_cap_host_h__
#ifdef HOST_CAPI_BUILD
#include "rtw_capi.h"
#include "rtw_modelmap.h"
#include "rtw_modelmap_simtarget.h"
typedef struct { rtwCAPI_ModelMappingInfo mmi ; } tryvm_host_DataMapInfo_T ;
#ifdef __cplusplus
extern "C" {
#endif
void tryvm_host_InitializeDataMapInfo ( tryvm_host_DataMapInfo_T * dataMap ,
const char * path ) ;
#ifdef __cplusplus
}
#endif
#endif
#endif
