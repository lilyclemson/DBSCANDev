
/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2019 HPCC Systems.  All rights reserved.
############################################################################## */IMPORT ML_Core;
IMPORT $.^ AS DBSCAN;
IMPORT $.datasets.frogDS_Small AS frog_data;

// Test to check the Num_clusters function

ds := frog_data.ds;

ML_Core.AppendSeqID(ds,id,dsID);
ML_Core.ToField(dsID,dsNF);

OUTPUT(DBSCAN.DBSCAN(0.3,10).Num_Clusters(dsNF));
OUTPUT(DBSCAN.DBSCAN(0.3,10).Num_Outliers(dsNF));