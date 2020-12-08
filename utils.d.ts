import SparseMap from 'mnemonist/sparse-map';
import {UndirectedLouvainIndex} from 'graphology-indices/neighborhood/louvain';

type PointerArray = Uint8Array | Uint16Array | Uint32Array | Float64Array;

type UndirectedLeidenAddendaOptions = {
  rng?: () => number;
};

export class UndirectedLeidenAddenda {
  B: number;
  C: number;
  resolution: number;
  randomness: number;
  index: UndirectedLouvainIndex;
  rng: () => number;

  communitiesOffsets: PointerArray;
  nodesSortedByCommunities: PointerArray;
  communitiesBounds: PointerArray;
  communityWeights: PointerArray;
  nonSingleton: PointerArray;
  externalEdgeWeightPerCommunity: PointerArray;
  belongings: PointerArray;
  neighboringCommunities: SparseMap<number>;
  cumulativeIncrement: Float64Array;

  constructor(index: UndirectedLeidenAddenda, options?: UndirectedLeidenAddendaOptions);
  groupByCommunities(): void;
  mergeNodesSubset(): void;
  refinePartition(): void;
}
