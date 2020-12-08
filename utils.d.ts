import SparseMap from 'mnemonist/sparse-map';
import {UndirectedLouvainIndex} from 'graphology-indices/neighborhood/louvain';

type PointerArray = Uint8Array | Uint16Array | Uint32Array | Float64Array;

export class UndirectedLeidenAddenda {
  B: number;
  resolution: number;
  index: UndirectedLouvainIndex;

  communitiesOffsets: PointerArray;
  nodesSortedByCommunities: PointerArray;
  communitiesBounds: PointerArray;
  clusterWeights: PointerArray;
  nonSingletonClusters: PointerArray;
  externalEdgeWeightPerCluster: PointerArray;
  belongings: PointerArray;
  neighboringCommunities: SparseMap<number>;
  cumulativeIncrement: Float64Array;

  groupByCommunities(): void;
  mergeNodesSubset(): void;
  refinePartition(): void;
}
