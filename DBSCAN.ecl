IMPORT ML_Core;
IMPORT ML_Core.Types AS Types;
IMPORT DBSCAN_Types AS Files;
IMPORT Std.system.Thorlib;


EXPORT DBSCAN(REAL8 eps = 0, UNSIGNED4 minPts = 2, STRING dist = 'euclidian', SET OF REAL8 dist_params = []):= MODULE

  EXPORT STREAMED DATASET(Files.l_stage3) locDBSCAN(STREAMED DATASET(Files.l_stage2) dsIn, //distributed data from stage 1
                                                    REAL8 eps = eps,   //distance threshold
                                                    UNSIGNED minPts = minPts, //the minimum number of points required to form a cluster,
                                                    STRING distance_func = dist,
                                                    SET OF REAL8 params = dist_params,
                                                    UNSIGNED4 localNode = Thorlib.node()
                                                    ) := EMBED(C++ : activity)

    #include<iostream>
    #include<bits/stdc++.h>

    using namespace std;

    string distanceFunc = "euclidian";
    vector<double> dist_params;

    struct node{
        uint16_t wi;
        uint64_t id;
        uint64_t nodeId;
        vector<float> fields;
        bool ifLocal = false;
        bool ifCore = false;
        node *parent = NULL;
        bool isModified = false;
        bool isVisited = false;
    };

    double euclidian(vector<float> a, vector<float> b){
        double sum=0;
        for(int i=0; i<a.size(); ++i){
            sum += (a[i] - b[i])*(a[i] - b[i]);
        }
        return sqrt(sum);
    }

    double haversine(vector<float> a, vector<float> b){
        int M=a.size();
        
        if(M!=2){
            return 0;
        }

        double lat1=a[0];
        double lat2=b[0];
        double lon1=a[1];
        double lon2=b[1];

        double sin_0 = sin(0.5 * (lat1 - lat2));
        double sin_1 = sin(0.5 * (lon1 - lon2));

        return (sin_0 * sin_0 + cos(lat1) * cos(lat2) * sin_1 * sin_1);
    }


    double manhattan(vector<float> a, vector<float> b){
        double ans=0;
        int M=a.size();
        
        for(int i=0;i<M;i++)
        ans=ans+(abs((a[i])-(b[i])));

        return ans;
    }

    double minkowski(vector<float> a, vector<float> b, int p){
        // sum(|x - y|^p)^(1/p)

        int m=a.size();
        double ans=0;
        for(int i=0;i<m;i++){
                ans+=pow(abs(a[i]-b[i]),p);
        }

        return pow(ans,(double)1/p);
    }

    double cosine(vector<float> a, vector<float> b){
        double ans=0, a1=0,a2=0;
        int m=a.size();
        for(int i=0;i<m;i++){
                ans+=a[i]*b[i];
                a1=a1+pow(a[i],2);
                a2=a2+pow(b[i],2);
        }
        return ans/(sqrt(a1)*sqrt(a2));
    }

    double chebyshev(vector<float> a, vector<float> b){
        // max(|x - y|)
        int m=a.size();
        float ans=0;
        for(int i=0;i<m;i++){
                ans=max(abs(a[i]-b[i]),ans);
        }
        return ans;
    }

    vector<node*> getNeighbors(vector<node*> ds, node *p, double eps){
        vector<node*> ret;
        for(uint64_t i=0; i<ds.size(); ++i){

            if(p->wi != ds[i]->wi) continue;

            double dist;

            if(distanceFunc.compare("cosine")==0)
                dist = cosine(p->fields, ds[i]->fields);
            else if(distanceFunc.compare("minkowski")==0)
                dist = minkowski(p->fields, ds[i]->fields, dist_params[0]);
            else if(distanceFunc.compare("manhattan")==0)
                dist = manhattan(p->fields, ds[i]->fields);
            else if(distanceFunc.compare("haversine")==0)
                dist = haversine(p->fields, ds[i]->fields);
            else if(distanceFunc.compare("chebyshev")==0)
                dist = chebyshev(p->fields, ds[i]->fields);
            else
                dist = euclidian(p->fields, ds[i]->fields);

            if(dist <= eps) ret.push_back(ds[i]);
        }
        return ret;
    }

    node* find(node *p){
        if(p->parent == NULL || p->parent == p){
            return p;
        } else {
            return find(p->parent);
        }
    }

    void Union(vector<node*> ds, node* a, node* b){
        node* x = find(a);
        node* y = find(b);
        if(x->ifCore && !y->ifCore){
            y->parent = x;
            return;
        } else if (!x->ifCore && y->ifCore){
            x->parent = y;
            return;
        }
        if(x == y) return;
        else if(x->id > y->id) y->parent = x;
        else x->parent = y;
    }

    void dbscan(vector<node*> ds, double eps, uint64_t minpts) {
        for(uint64_t i=0; i<ds.size(); ++i){

            if(!ds[i]->ifLocal) continue;

            ds[i]->isModified = true;
            vector<node*> neighs = getNeighbors(ds,ds[i],eps);

            if(neighs.size() >= minpts){
                ds[i]->ifCore = true;
                ds[i]->parent = NULL;
                for(uint64_t n=0; n < neighs.size(); ++n){
                    neighs[n]->isModified = true;
                    if(neighs[n]->ifLocal){
                        if(neighs[n]->ifCore){
                            Union(ds,ds[i],neighs[n]);
                        } else {
                            if(neighs[n]->isVisited) continue;
                            neighs[n]->isVisited = true;
                            Union(ds,ds[i],neighs[n]);
                        }
                    } else {
                        if(!neighs[n]->ifCore){
                            if(getNeighbors(ds,neighs[n],eps).size()>=minpts){
                                neighs[n]->ifCore = true;
                                Union(ds,ds[i],neighs[n]);
                            } else {
                                if(neighs[n]->isVisited) continue;
                                neighs[n]->isVisited = true;
                                Union(ds,ds[i],neighs[n]);
                            }
                        } else {
                            Union(ds,ds[i],neighs[n]);
                        }
                    }
                }
            }
        }
    }

    struct retRecord{
    uint32_t wi;
    uint32_t id;
    uint32_t parentId;
    uint32_t nodeId;
    bool if_local;
    bool if_core;
    };

    class ResultStream : public RtlCInterface, implements IRowStream {
    public:
        ResultStream(IEngineRowAllocator *_ra, IRowStream *_ds, uint64_t minpts, double eps, unsigned long long lnode) : ra(_ra), ds(_ds), Lnode(lnode){
            byte* p;
            count = 0;
            while((p=(byte*)ds->nextRow())){
                node temp;
                temp.wi = *((uint16_t*)p); p += sizeof(uint16_t);
                temp.id = *((unsigned long long*)p); p += 2*sizeof(unsigned long long);
                temp.nodeId = *((unsigned long long*)p); p += sizeof(unsigned long long);
                p += sizeof(bool);
                uint32_t lenFields = *((uint32_t*)p); 
                p += sizeof(uint32_t);
                for(uint32_t i=0; i<lenFields/sizeof(float); ++i){
                    float f = *((float*)p); p += sizeof(float);
                    temp.fields.push_back(f);
                }
                temp.ifLocal = (lnode == temp.nodeId); p += sizeof(bool);
                temp.ifCore = false;
                dataset.push_back(temp);
            }

            for(uint64_t i=0; i<dataset.size(); ++i){
                results.push_back(&dataset[i]);
            }

            dbscan(results, eps, minpts);
        }

        RTLIMPLEMENT_IINTERFACE
        virtual const void* nextRow() override {
            RtlDynamicRowBuilder rowBuilder(ra);
            if(count < results.size()){
                uint32_t lenRec = 4*sizeof(uint32_t) + 2*sizeof(bool);
                byte* p = (byte*)rowBuilder.ensureCapacity(lenRec, NULL);

                while(!dataset[count].isModified && count < dataset.size())
                    count++;
                
                if(count >= results.size()) return NULL;
                
                *((uint32_t*)p) = dataset[count].wi; p += sizeof(uint32_t);
                *((uint32_t*)p) = Lnode; p += sizeof(uint32_t);
                *((uint32_t*)p) = dataset[count].id; p += sizeof(uint32_t);
                *((uint32_t*)p) = find(&dataset[count])->id; p += sizeof(uint32_t);
                *((bool*)p) = dataset[count].ifLocal; p += sizeof(bool);
                *((bool*)p) = dataset[count].ifCore;

                count++;
                return rowBuilder.finalizeRowClear(lenRec);
            } else {
                return NULL;
            }
        }

        virtual void stop() override{}

    protected:
        Linked<IEngineRowAllocator> ra;
        unsigned count;
        vector<node> dataset;
        vector<node*> results;
        IRowStream *ds;
        uint64_t Lnode;
    };

    #body

    distanceFunc = distance_func;
    double* p = (double*)params;

    for(uint32_t i=0; i<lenParams/sizeof(double); ++i){
        dist_params.push_back(*p);
        p += sizeof(double);
    }

    transform(distanceFunc.begin(),distanceFunc.end(),distanceFunc.begin(),::tolower);

    return new ResultStream(_resultAllocator, dsin, minpts, eps, localnode);
  ENDEMBED;//end locDBSCAN

  //Layout for Ultimate() and Loop_Func()
  EXPORT l_ultimate := RECORD
    UNSIGNED4 wi;
    UNSIGNED4 id;
    UNSIGNED4 parentID;
    UNSIGNED4 largestID := 0;
    UNSIGNED4 ultimateID := 0;
  END;

  EXPORT STREAMED DATASET(l_ultimate) ultimate(STREAMED DATASET(l_ultimate) dsin, UNSIGNED4 pointcount) := EMBED(C++:activity)
    #include <stdio.h>
    struct upt
    {
      uint32_t wi;
      uint32_t id;
      uint32_t pid;
    };

    class MyStreamInlineDataset : public RtlCInterface, implements IRowStream
      {
        public:
            MyStreamInlineDataset(IEngineRowAllocator * _resultAllocator, IRowStream * _ds, uint32_t _pc)
              :resultAllocator(_resultAllocator), ds(_ds), pc(_pc)
              {
                uptable = (upt*) rtlMalloc(pc * sizeof(upt));
                for(uint32_t i = 0; i < pc; i++)
                {
                  uptable[i].wi = 0;
                  uptable[i].id = 0;
                  uptable[i].pid = 0;
                };
                calculated = false;
                rc = 0;
                lastgroupend = 0;
                curWi = 0;
                lastid = 0;
              }
              ~MyStreamInlineDataset(){
              // rtlFree(uptable);
              }

            RTLIMPLEMENT_IINTERFACE
            //calculate the ultimate id
            virtual const void *nextRow() override
            {
                if(!calculated){
                    while(true)
                    {
                        const byte * next = (const byte *)ds->nextRow();
                        if (!next) break;
                        const byte * pos = next;
                        uint32_t wi = *(uint32_t*)pos;
                        pos += sizeof(uint32_t);
                        uint32_t id = *(uint32_t*)pos;
                        pos += sizeof(uint32_t);
                        uint32_t pid = *(uint32_t *) pos;
                        if(curWi == 0){
                          curWi = wi;
                        }
                        if(curWi != wi){
                          curWi = wi;
                          lastgroupend = lastid;
                        }
                        id += lastgroupend;
                        pid += lastgroupend;
                        if (id > 0 && id <= pc)
                        {
                        uptable[id -1].wi = wi;
                        uptable[id -1].id = id;
                        uptable[id -1].pid = pid;
                        }
                        lastid = id;
                        rtlReleaseRow(next);
                    }// End while()

                    for(uint32_t i = 0; i < pc; i++)
                    {
                      uint32_t wi = uptable[i].wi;
                      uint32_t id = uptable[i].id;
                      uint32_t pid = uptable[i].pid;
                      if(id == 0) continue;
                      while(id != pid)
                      {
                        id = pid;
                        if(pid - 1  >= pc){
                          break;
                        }
                        if(uptable[pid -1].pid == 0 || uptable[pid -1].wi != wi){
                          break;
                        }else{
                          pid = uptable[pid -1].pid;
                        }
                      }
                      uptable[i].pid = pid;
                    };// end for()

                  calculated = true;
                  lastgroupend = 0;
                  curWi = 0;
                  lastid = 0;
                }//end if(!calculated)

                byte* row;
                RtlDynamicRowBuilder rowBuilder(resultAllocator);
                uint32_t returnsize = 5*sizeof(uint32_t);
                while(rc < pc && uptable[rc].id == 0){ rc++;}
                if(rc < pc)
                {
                  row = rowBuilder.ensureCapacity(returnsize, NULL);
                  void * pos = row;
                  uint32_t id = uptable[rc].id;
                  uint32_t pid = uptable[rc].pid;
                  uint32_t wi = uptable[rc].wi;

                  if(curWi == 0){
                    curWi = wi;
                  }
                  if(curWi != wi){
                    curWi = wi;
                    lastgroupend = lastid;
                  }
                  id = id - lastgroupend;
                  pid = pid - lastgroupend;
                  *(uint32_t *)pos = wi;
                  pos += sizeof(uint32_t);
                  *(uint32_t *)pos = id;
                  pos += sizeof(uint32_t);
                  *(uint32_t *)pos = pid;
                  pos += sizeof(uint32_t);
                  *(uint32_t *)pos = rc;
                  pos += sizeof(uint32_t);
                  *(uint32_t *)pos = lastgroupend;
                  lastid = id;
                  rc++;
                  return rowBuilder.finalizeRowClear(returnsize);
                }else{
                  return NULL;
                }// end if()

            }// end nextRow()

            virtual void stop() override
            {
                // ds->stop();
            }

            protected:
                Linked<IEngineRowAllocator> resultAllocator;
                IRowStream * ds;
                uint32_t pc;
                upt * uptable;
                bool calculated;
                uint32_t rc;// row counter
                uint32_t lastgroupend;
                uint32_t curWi;
                uint32_t lastid;
      };

    #body

          return new MyStreamInlineDataset(_resultAllocator, dsin, pointcount);

  ENDEMBED;//end ultimate()

  //LOOP to get the final result/ultimateID
  EXPORT Loop_Func(DATASET(l_ultimate) ds, UNSIGNED c) := FUNCTION
        tempLayout := RECORD
          UNSIGNED4 wi;
          UNSIGNED4 id;
          UNSIGNED4 newParentID;
        END;
        tempChanges := PROJECT(ds, TRANSFORM(tempLayout,
                                      SELF.wi := LEFT.wi,
                                      SELF.id := LEFT.ultimateID,
                                      SELF.newParentID := LEFT.largestID), LOCAL);
        changes := DEDUP(SORT(tempChanges, wi, id, -newParentID, LOCAL), wi, id, LOCAL);
        newParent := JOIN(ds, changes, LEFT.wi = RIGHT.wi AND LEFT.id = RIGHT.id, TRANSFORM(RECORDOF(LEFT),
                                                              SELF.parentID := IF(right.id > 0, RIGHT.newParentID, LEFT.parentID),
                                                              SELF := LEFT), LEFT OUTER, LOCAL);

        newUltimate :=  Ultimate(newParent, c);
        rst := JOIN(newParent, newUltimate, LEFT.wi = RIGHT.wi AND LEFT.id = RIGHT.id, TRANSFORM(l_ultimate,
                                                          SELF.ultimateID := RIGHT.parentID,
                                                          SELF := LEFT));
        RETURN rst;
  END;//end loop_func()

  /**
    * fit() function i
    */
  EXPORT DATASET(ML_Core.Types.ClusterLabels) fit(DATASET(Types.NumericField) ds) := FUNCTION

    /**
      * Stage 1: Transform and distribute input dataset ds for local clustering in stage 2.
      */

    //Evenly distribute the data
    Xnf1 := DISTRIBUTE(ds, id);
    //Transform to 1_stage1
    X0 := PROJECT(Xnf1, TRANSFORM(
                                  Files.l_stage1,
                                  SELF.fields := [LEFT.value],
                                  SELF.nodeId := Thorlib.node(),
                                  SELF := LEFT),
                                  LOCAL);
    X1 := SORT(X0, wi, id, number, LOCAL);
    X2 := ROLLUP(X1, TRANSFORM(
                                Files.l_stage1,
                                SELF.fields := LEFT.fields + RIGHT.fields,
                                SELF := LEFT),
                                wi, id,
                                LOCAL);
    //Transform to l_stage2
    X3 := PROJECT(X2, TRANSFORM(
                                Files.l_stage2,
                                SELF.parentID := LEFT.id,
                                SELF := LEFT),
                                LOCAL);
    //Braodcast for local clustering.
    X := DISTRIBUTE(X3, ALL);

    /**
      * Stage 2: local clustering on each node
      */

    rds := locDBSCAN(X);

    /**
      * Stage 3: global merge the local clustering results to the final clustering result
      */


    //get non_outliers and its largest parentID
    rds1 := rds( NOT( if_core = FALSE AND id = parentID));
    non_outliers := DEDUP(SORT(rds1,wi, id,-parentID),wi,id );

    //Get outliers
    outliers := PROJECT(JOIN(rDS, non_outliers, LEFT.wi = RIGHT.wi AND LEFT.id = RIGHT.id, LEFT ONLY), TRANSFORM(l_ultimate,
                                                                                  SELF.ultimateid := LEFT.parentid,
                                                                                  SELF := LEFT));
    unfiltered := rDS(if_local = TRUE);
    ntunfiltered := COUNT(unfiltered );
    dds := DISTRIBUTE(unfiltered, wi); //
    f0 := PROJECT(NOCOMBINE(dds), TRANSFORM({l_Ultimate, UNSIGNED4 nodeid}, SELF.nodeid := Thorlib.node(), SELF := LEFT));
    t := TABLE(f0, { nodeid , cnt := COUNT(GROUP)}, nodeid, LOCAL);
    c := t(nodeid = thorlib.node())[1].cnt;

    //get local core points
    f1 := rDS(if_local = TRUE AND if_core=TRUE);
    f2 := DISTRIBUTE(f1, wi);
    localCores := SORT(PROJECT(NOCOMBINE(f2),TRANSFORM(l_ultimate, SELF := LEFT), LOCAL), wi, id, LOCAL);

    locals_ultimate:=  ultimate(localCores, c);// all the ultimates for locals

    //get largestID for the core points
    largest := DISTRIBUTE(non_outliers(if_core = TRUE), wi);

    //Prepare the input dataset 'initial' for Loop_Func()
    //Join largest and locals_ultimate
    initial0 := JOIN(largest, locals_ultimate,
                    LEFT.wi = RIGHT.wi
                    AND
                    LEFT.id = RIGHT.id,
                    TRANSFORM(l_ultimate,
                              SELF.ultimateID := RIGHT.parentID,
                              SELF.largestID := LEFT.parentID,
                              SELF := LEFT), LOCAL);

    //Join locals
    initial := JOIN(initial0, localCores,
                    LEFT.wi = RIGHT.wi
                    AND
                    LEFT.id = RIGHT.id,
                    TRANSFORM(l_ultimate,
                            SELF.parentID := RIGHT.parentID,
                            SELF := LEFT), LOCAL);

    l := LOOP(initial, LEFT.id > 0, EXISTS(ROWS(LEFT)(ultimateID < largestID)), LOOP_Func(ROWS(LEFT), c) );
    //Update the parentID of all non_outliers from the result
    update_non_outliers := JOIN(non_outliers, l, LEFT.wi = RIGHT.wi AND LEFT.parentid = RIGHT.id, TRANSFORM(l_ultimate,
                                                                    SELF.ultimateID := IF(right.id =0, LEFT.parentid, RIGHT.ultimateID),
                                                                    SELF:= LEFT),
                                                                    LEFT OUTER);
    //Convert ultimate id of outliers to zero
    outliers1 := PROJECT(outliers, TRANSFORM(ML_Core.Types.ClusterLabels,
                                             SELF.label := 0,
                                             SELF.id := LEFT.id,
                                             SELF.wi := LEFT.wi));
    //Get mapping to 1-based cluster labels
    map0 := TABLE(update_non_outliers,{wi,ultimateId,id_min:=MIN(GROUP,id)},wi,ultimateId);
    map1 := PROJECT(SORT(map0,wi,id_min), TRANSFORM(RECORDOF(map0),
                                    SELF.wi:=LEFT.wi,
                                    SELF.ultimateId:=LEFT.ultimateID,
                                    SELF.id_min:=COUNTER));
    //Map 1 based indices to the cluster labels
    result0 := PROJECT(update_non_outliers, TRANSFORM(ML_Core.Types.ClusterLabels,
                                                      SELF.label := map1(wi=LEFT.wi and ultimateid=LEFT.ultimateId)[1].id_min
                                                                    - MIN(map1(wi=LEFT.wi),map1(wi=LEFT.wi).id_min)
                                                                    + 1,
                                                      SELF.wi := LEFT.wi,
                                                      SELF.id := LEFT.id));
    //Combine outliers and non outliers to form final result
    result1 := result0 + outliers1;
    RETURN result1;
  END;//end fit()
  /**
    * Num_Clusters
    *
    * Provides the number of clusters that the given dataset will be divided into
    * when clustered by the DBSCAN algorithm.
    *
    * @param ds The dataset to be used for the clustering operation in DATASET(NumericField) format
    * @return DATASET(l_num_clusters) The number of clusters, per work item.
    * @see DBSCAN_Types.l_num_clusters
    */
  EXPORT DATASET(Files.l_num_clusters) Num_Clusters(DATASET(Types.NumericField) ds) := FUNCTION
    //Find clustering of ds
    clustering := Fit(ds);
    //Find maxmimum label of X samples per work item
    result0 := TABLE(clustering,{wi,num:=MAX(GROUP,label)},wi);
    //Project to match return type
    result1 := PROJECT(result0, TRANSFORM(Files.l_num_clusters,
                                          SELF.wi := LEFT.wi,
                                          SELF.num := LEFT.num));
    RETURN result1;
  END;//end Num_Clusters()
  /**
    * Num_Outliers
    *
    * Provides the number of outliers that the given dataset will have
    * when clustered by the DBSCAN algorithm.
    *
    * @param ds The dataset to be used for the clustering operation in DATASET(NumericField) format
    * @return DATASET(l_num_clusters) The number of outliers, per work item.
    * @see DBSCAN_Types.l_num_clusters
    */
  EXPORT DATASET(Files.l_num_clusters) Num_Outliers(DATASET(Types.NumericField) ds) := FUNCTION
    //Find clustering of ds
    clustering := Fit(ds);
    //Find number of outliers per work item
    outliers := TABLE(clustering(label=0),{wi,num:=COUNT(GROUP)},wi);
    //Project to match return type
    result := PROJECT(outliers, TRANSFORM(Files.l_num_clusters,
                                          SELF.wi := LEFT.wi,
                                          SELF.num := LEFT.num));
    RETURN result;
  END;//end Num_Outliers()

END;//end DBSCAN
