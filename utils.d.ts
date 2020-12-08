type PointerArray = Uint8Array | Uint16Array | Uint32Array | Float64Array;

export class UndirectedLeidenAddenda {
  B: number;

  communitiesOffsets: PointerArray;
  nodesSortedByCommunities: PointerArray;
  communitiesBounds: PointerArray;
}
