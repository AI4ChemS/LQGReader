from __future__ import annotations

import numpy as np

from LQGReader.json_reader import Tiling


class QuotientVertex:
    def __init__(self, rank: int, position=None, attributes=None) -> NoneType:
        self.rank = rank
        self.position = position
        self.attributes = attributes

        self.links


class QuotientLink:
    def __init__(self, start_vert, start_shift, end_vert, end_shift, join_vert=None, join_shift=None) -> NoneType:
        self.start_vert = start_vert
        self.start_shift = start_shift

        self.end_vert = end_vert
        self.end_shift = end_shift

        self.join_vert = join_vert
        self.join_shift = join_shift


class QuotientGraph:
    def __init__(self, vertices: list[QuotientVertex], links: dict[str, list[QuotientLink]]) -> NoneType:
        self.vertices = vertices
        self.links = links
    
    @staticmethod
    def build_from_tiling(tiling: Tiling) -> QuotientGraph:
        vertices = []
        for node in tiling.node_list.node_list:
            new_vertex = QuotientVertex(
                rank=1,
                position=node.get_position(),
            )
            




