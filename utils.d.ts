import SparseMap from 'mnemonist/sparse-map';
import {UndirectedLouvainIndex} from 'graphology-indices/neighborhood/louvain';

type PointerArray = Uint8Array | Uint16Array | Uint32Array | Float64Array;

type UndirectedLeidenAddendaOptions = {
  rng?: () => number;
};

export class UndirectedLeidenAddenda {
  B: number;
  resolution: number;
  randomness: number;
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

  constructor(index: UndirectedLeidenAddenda, options?: UndirectedLeidenAddendaOptions);
  groupByCommunities(): void;
  mergeNodesSubset(): void;
  refinePartition(): void;
}