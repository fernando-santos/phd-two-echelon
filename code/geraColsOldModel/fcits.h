/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_FCITS
#define _H_FCITS

typedef struct
{
  int A;
  int B;
  double X;
} EdgeRec;
typedef EdgeRec *EdgePtr;

typedef struct
{
  int PartitionSize;
  ReachPtr PartitionPtr;
  int EPtrSize;
  int EPtrDim;
  EdgePtr EPtr;
  int ShrunkEdgeNr;
} TreeSearchRec;
typedef TreeSearchRec *TreeSearchPtr;

void FCITS_ComputeLHS(ReachPtr CompsRPtr,
                      int NoOfComps,
                      int NoOfSuperNodes,
                      double **FlowMatrix,
                      ReachPtr FlowRPtr,
                      double *LHS);

void FCITS_ComputeRHS(ReachPtr PartitionPtr,
                      int PartitionSize,
                      int *SuperNodeDemand,
                      int CAP,
                      double *RHS);

void FCITS_ShrinkPartition(int NoOfSuperNodes,
                           ReachPtr FlowPtr,
                           double **FlowMatrix,
                           ReachPtr PartitionPtr,
                           int PartitionSize,
                           int NodeA,
                           int NodeB,
                           ReachPtr NewPartitionPtr,
                           int *NewPartitionSize);

void FCITS_CreateEPtrForPartition(int NoOfSuperNodes,
                                  ReachPtr FlowPtr,
                                  double **FlowMatrix,
                                  ReachPtr PartitionPtr,
                                  int PartitionSize,
                                  int Level,
                                  TreeSearchPtr TreePtr);

void FCITS_CheckForDominance(int NoOfSuperNodes,
                             int CurrentLevel,
                             TreeSearchPtr TreePtr,
                             char *Dominated);

void FCITS_TreeSearch(int NoOfSuperNodes,
                      int CAP,
                      int *SuperNodeDemand,
                      double **FlowMatrix,
                      ReachPtr FlowRPtr,
                      int MaxCuts,
                      int MaxFCITSLoops,
                      int *GeneratedCuts,
                      double *MaxViolation,
                      double *CutsRHS,
                      ReachPtr CutsRPtr);

void FCITS_ComputeFlowMatrix(ReachPtr SupportPtr,
                             int NoOfCustomers,
                             double **XMatrix,
                             ReachPtr SuperNodesRPtr,
                             int NoOfSuperNodes,
                             double **FlowMatrix);

void FCITS_MainCutGen(ReachPtr SupportPtr,
                      int NoOfCustomers,
                      int *Demand,
                      int CAP,
                      double **XMatrix,
                      ReachPtr InitSuperNodesRPtr,
                      ReachPtr InitSAdjRPtr,
                      int *InitSuperDemand,
                      int InitShrunkGraphCustNodes,
                      int MaxFCITSLoops,
                      int MaxNoOfCuts,
                      double *MaxViolation,
                      int *NoOfGeneratedCuts,
                      CnstrMgrPointer CutsCMP);

#endif

