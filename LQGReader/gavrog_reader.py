"""
Positions are always relative, unless specified to be global. 
"""
from __future__ import annotations

import itertools
import json

import numpy as np

from LQGReader import utils


class Basis:
    def __init__(self, vectors: np.ndarray) -> NoneType:
        self.vectors = vectors
        self.op = self.vectors.T
        self.inv_op = np.linalg.inv(self.op)
    
    def rel_to_global(self, positions: np.ndarray) -> np.ndarray:
        global_positions = np.dot(self.op, positions.T).T

        return global_positions
    
    def global_to_rel(self, global_positions: np.ndarray) -> np.ndarray:
        rel_positions = np.dot(self.inv_op, global_positions.T).T

        return rel_positions
    
    def global_norm(self, x: np.ndarray):
        global_x = self.rel_to_global(x)
        global_norm = np.linalg.norm(global_x)

        return global_norm
    

class Node:
    def __init__(self, unit_cell_position: np.ndarray, node_id=None) -> NoneType:
        self.unit_cell_position = unit_cell_position
        self.id = node_id

        self.shift_nodes = []
    
    @staticmethod 
    def build_from_position(position: np.ndarray, node_id=None) -> Node:
        unit_cell_position = utils.compute_unit_cell_position(position)

        new_node = Node(
            unit_cell_position=unit_cell_position,
            node_id=node_id
        )

        return new_node
    
    def equivalent(self, node: Node, basis: Basis, epsilon=1e-2) -> bool:
        periodic_diff = utils.periodic_diff(
            self.unit_cell_position,
            node.unit_cell_position, 
        )
        global_norm = basis.global_norm(periodic_diff)

        return global_norm < epsilon
    
    def __repr__(self,) -> str:
        str_rep = "Node. "
        if self.id is not None:
            str_rep += f"ID: {self.id}, "
        
        str_rep += f"Position: {np.round(self.unit_cell_position, decimals=3)}"

        return str_rep
    
    def get_position(self,) -> np.ndarray:
        return self.unit_cell_position


class NodeList:
    def __init__(self,) -> NoneType:
        self.node_list = []
    
    def find(self, find_node: Node, basis: Basis, epsilon=1e-2) -> Node:
        for node in self.node_list:
            if find_node.equivalent(node, basis=basis, epsilon=epsilon):
                return node

        return None

    def add_from_position(self, position: np.ndarray, basis: Basis, epsilon=1e-2) -> NoneType:
        new_node = Node.build_from_position(
            position=position,
            node_id=len(self.node_list),
        )
        found_node = self.find(
            find_node=new_node,
            basis=basis,
            epsilon=epsilon,
        )

        if found_node is not None:
            return
        
        self.node_list.append(new_node)
    
    def add_from_multiple_positions(self, positions: np.ndarray, basis: Basis, epsilon=1e-2) -> NoneType:
        for pos in positions:
            self.add_from_position(position=pos, basis=basis, epsilon=epsilon)
    
    @staticmethod
    def build_from_global_positions(global_positions: np.ndarray, basis: Basis, epsilon=1e-2) -> NodeList:
        positions = basis.global_to_rel(global_positions)

        print(global_positions)
        print(positions)

        new_node_list = NodeList()
        new_node_list.add_from_multiple_positions(
            positions=positions,
            basis=basis,
            epsilon=epsilon,
        )

        return new_node_list
    
    @staticmethod
    def build_from_json(face_dict_list: dict, basis: Basis, epsilon=1e-2):
        global_positions = list(
            itertools.chain.from_iterable(
                [face["nodePositions"] for face in face_dict_list]
            )
        )
        global_positions = np.array(global_positions)

        node_list = NodeList.build_from_global_positions(
            global_positions=global_positions,
            basis=basis,
            epsilon=epsilon,
        )

        return node_list
    
    def __repr__(self,) -> str:
        str_rep = "Node List\n"
        for node in self.node_list:
            str_rep += str(node) + '\n'
        
        return str_rep.strip()


class ShiftNode:
    def __init__(self, node: Node, shift: np.ndarray) -> NoneType:
        self.node = node
        self.shift = shift

        self.shift_edges = []

        self.node.shift_nodes.append(self)

    @staticmethod
    def build_from_position(position: np.ndarray, node_list: NodeList, basis: Basis, epsilon=1e-2) -> ShiftNode:
        temp_node = Node.build_from_position(position=position)

        node = node_list.find(
            find_node=temp_node,
            basis=basis,
            epsilon=epsilon
        )

        if node is None:
            raise ValueError(f"No node found for: {temp_node}")

        _, shift = utils.compute_shift_position(position=position)
        new_shift_node = ShiftNode(
            node=node,
            shift=shift,
        )

        return new_shift_node
    
    def get_position(self,) -> np.ndarray:
        return self.node.get_position() + self.shift
    
    def __repr__(self,) -> str:
        str_rep = f"Shift Node. Shift: {self.shift}. {self.node}"

        return str_rep
    
    def __eq__(self, shift_node: ShiftNode):
        return self.node == shift_node.node


class Edge:
    def __init__(self, shift_node_set: set[ShiftNode], edge_id=None) -> NoneType:
        self.shift_node_set = shift_node_set
        self.id = edge_id

        self.shift_edges = []
    
    @staticmethod
    def build_from_positions(positions: np.ndarray, node_list: NodeList, basis: Basis, epsilon=1e-2, edge_id=None) -> Edge:
        shift_node_set = set(
            [ShiftNode.build_from_position(pos, node_list=node_list, basis=basis, epsilon=epsilon) for pos in positions]
        )
        
        new_edge = Edge(
            shift_node_set={shift_node_A, shift_node_B},
            edge_id=edge_id,
        )

        return new_edge
    
    def equivalent(self, edge: Edge, basis: Basis, epsilon=1e-2):

    
    def get_position(self,) -> np.ndarray:
        edge_center = np.mean(
            np.array(
                [shift_node.get_position() for shift_node in self.shift_node_set]
            ),
            axis=0,
        )

        return edge_center
    
    def __repr__(self,) -> str:
        str_rep = f"Edge. "
        if self.id is not None:
            str_rep += f"ID: {self.id}, "
        
        str_rep += f"Center: {np.round(self.edge_center.unit_cell_position, decimals=3)}, {self.vector_list}"

        return str_rep


class EdgeList:
    def __init__(self,) -> NoneType:
        self.edge_list = []
    
    def find(self, find_edge: Edge, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> Edge:
        for edge in self.edge_list:
            if find_edge.equivalent(edge, basis, epsilon_center=epsilon_center, epsilon_vector=epsilon_center):
                return edge
        
        return None
    
    def add_from_positions(self, pos_A: np.ndarray, pos_B: np.ndarray, node_list: NodeList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> NoneType:
        new_edge = Edge.build_from_positions(
            pos_A=pos_A,
            pos_B=pos_B,
            node_list=node_list,
            basis=basis,
            edge_id=len(self.edge_list),
        )
        found_edge = self.find(
            find_edge=new_edge,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector
        )

        if found_edge is not None:
            return
        
        self.edge_list.append(new_edge)
    
    def add_from_multiple_positions(self, edge_positions: np.ndarray, node_list: NodeList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> NoneType:
        for (pos_A, pos_B) in edge_positions:
            self.add_from_positions(
                pos_A=pos_A,
                pos_B=pos_B,
                node_list=node_list,
                basis=basis,
                epsilon_center=epsilon_center,
                epsilon_vector=epsilon_vector
            )
    
    @staticmethod
    def build_from_global_positions(global_edge_positions: np.ndarray, node_list: NodeList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> EdgeList:
        edge_positions = basis.global_to_rel(global_positions=global_edge_positions)

        formatted_edge_positions = []
        for i in range(0, edge_positions.shape[0], 2):
            formatted_edge_positions.append(
                [edge_positions[i], edge_positions[i + 1]]
            )
        formatted_edge_positions = np.array(formatted_edge_positions)

        new_edge_list = EdgeList()
        new_edge_list.add_from_multiple_positions(
            edge_positions=formatted_edge_positions,
            node_list=node_list,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
        )

        return new_edge_list
    
    @staticmethod
    def build_from_json(face_dict_list: dict, basis: Basis, node_list: NodeList, epsilon_center=1e-2, epsilon_vector=1e-2,) -> EdgeList:
        global_position_list = []
        for face in face_dict_list:
            global_node_positions = face["nodePositions"]
            num_nodes = len(global_node_positions)
            for index in range(num_nodes):
                global_position_list.append(global_node_positions[index])
                global_position_list.append(global_node_positions[(index + 1) % num_nodes])
        
        global_positions = np.array(global_position_list)

        edge_list = EdgeList.build_from_global_positions(
            global_edge_positions=global_positions,
            node_list=node_list,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
        )

        return edge_list
    
    def __repr__(self,) -> str:
        str_rep = "Edge List:\n"
        for edge in self.edge_list:
            str_rep += str(edge) + '\n'
        
        return str_rep.strip()


class ShiftEdge:
    def __init__(self, edge: Edge, shift: np.ndarray) -> NoneType:
        self.edge = edge
        self.shift = shift

        self.shift_faces = []

        self.edge.shift_edges.append(self)
        self.edge.shift_node_A.shift_edges.append(self)
        self.edge.shift_node_B.shift_edges.append(self)

    @staticmethod
    def build_from_positions(pos_A: np.ndarray, pos_B, node_list: NodeList, edge_list: EdgeList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> ShiftEdge:
        temp_edge = Edge.build_from_positions(
            pos_A=pos_A,
            pos_B=pos_B,
            node_list=node_list,
            basis=basis,
        )

        edge = edge_list.find(
            find_edge=temp_edge,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
        )

        if edge is None:
            raise ValueError(f"No edge found for: {temp_edge}")

        _, shift = utils.compute_shift_position(position=(pos_A + pos_B) / 2.0)
        new_shift_edge = ShiftEdge(
            edge=edge,
            shift=shift,
        )

        return new_shift_edge
    
    def get_position(self,) -> np.ndarray:
        return self.edge.get_position() + self.shift
    
    def __repr__(self,) -> str:
        str_rep = f"Shift Node. Shift: {self.shift}. {self.edge}"

        return str_rep


class Face:
    def __init__(self, face_center: Node, vector_list: VectorList, shift_edge_list: list[ShiftEdge], face_id=None) -> NoneType:
        self.face_center = face_center
        self.vector_list = vector_list
        self.shift_edge_list = shift_edge_list
        self.id = face_id

        self.shift_faces = []
    
    @staticmethod
    def build_from_positions(positions: np.ndarray, node_list: NodeList, edge_list: EdgeList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2, face_id=None):
        face_center_pos, center_shift = utils.compute_shift_position(
            np.mean(positions, axis=0)
        )
        face_center = Node(face_center_pos)
        
        shifted_positions = positions - center_shift
        vector_list = VectorList()
        for shifted_pos in shifted_positions:
            vector_list.add(shifted_pos - face_center.unit_cell_position)

        shift_edge_list = []
        num_nodes = positions.shape[0]
        for index in range(num_nodes):
            shift_edge_list.append(
                ShiftEdge.build_from_positions(
                    pos_A=positions[index],
                    pos_B=positions[(index + 1) % num_nodes],
                    node_list=node_list,
                    edge_list=edge_list,
                    basis=basis,
                    epsilon_center=epsilon_center,
                    epsilon_vector=epsilon_vector
                )
            )
        
        new_face = Face(
            face_center=face_center,
            vector_list=vector_list,
            shift_edge_list=shift_edge_list,
            face_id=face_id,
        )

        return new_face

    def equivalent(self, face: Face, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> bool:
        equiv_centers = self.face_center.equivalent(
            face.face_center, 
            basis=basis, 
            epsilon=epsilon_center,
        )

        equiv_vectors = self.vector_list.equivalent(
            face.vector_list, 
            basis=basis, 
            epsilon=epsilon_vector,
        )

        return equiv_centers and equiv_vectors
    
    def get_position(self,) -> np.ndarray:
        return self.face_center.get_position()

    def __repr__(self,) -> str:
        str_rep = "Face. "
        if self.id is not None:
            str_rep += f"ID: {self.id}, "
        
        str_rep += f"Center: {np.round(self.face_center.unit_cell_position, decimals=3)}, {self.vector_list}"
        
        return str_rep


class FaceList:
    def __init__(self,) -> NoneType:
        self.face_list = []
    
    def find(self, find_face: Face, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> Face:
        for face in self.face_list:
            if find_face.equivalent(face, basis, epsilon_center=epsilon_center, epsilon_vector=epsilon_vector):
                return face
        
        return None
    
    def add_from_positions(self, positions: np.ndarray, node_list: NodeList, edge_list: EdgeList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> NoneType:
        new_face = Face.build_from_positions(
            positions=positions,
            node_list=node_list,
            edge_list=edge_list,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
            face_id=len(self.face_list)
        )
        found_face = self.find(
            find_face=new_face,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
        )

        if found_face is not None:
            return
        
        self.face_list.append(new_face)
    
    @staticmethod
    def build_from_json(face_dict_list: dict, basis: Basis, node_list: NodeList, edge_list: EdgeList, epsilon_center=1e-2, epsilon_vector=1e-2,) -> FaceList:
        face_list = FaceList()

        for face in face_dict_list:
            global_node_positions = np.array(face["nodePositions"])
            positions = basis.global_to_rel(global_node_positions)

            face_list.add_from_positions(
                positions=positions,
                node_list=node_list,
                edge_list=edge_list,
                basis=basis,
                epsilon_center=epsilon_center,
                epsilon_vector=epsilon_vector,
            )

        return face_list
    
    def __repr__(self,) -> str:
        str_rep = "FaceList\n"
        for face in self.face_list:
            str_rep += str(face) + '\n'

        return str_rep


class ShiftFace:
    def __init__(self, face, shift) -> NoneType:
        self.face = face
        self.shift = shift

        self.tiles = []

        self.face.shift_faces.append(self)
        for index in range(len(self.face.shift_edge_list)):
            self.face.shift_edge_list[index].shift_faces.append(self)

    @staticmethod
    def build_from_positions(positions: np.ndarray, node_list: NodeList, edge_list: EdgeList, face_list: FaceList, basis: Basis, epsilon_center=1e-2, epsilon_vector=1e-2) -> ShiftFace:
        temp_face = Face.build_from_positions(
            positions=positions,
            node_list=node_list,
            edge_list=edge_list,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
        )

        face = face_list.find(
            find_face=temp_face,
            basis=basis,
            epsilon_center=epsilon_center,
            epsilon_vector=epsilon_vector,
        )

        if face is None:
            raise ValueError(f"No face found for: {temp_face}")
        
        _, shift = utils.compute_shift_position(np.mean(positions, axis=0))
        new_shift_face = ShiftFace(
            face=face,
            shift=shift,
        )

        return new_shift_face
    
    def get_position(self,) -> np.ndarray:
        return self.face.get_position() + self.shift
    
    def __repr__(self,) -> str:
        str_rep = f"Shift Face. Shift: {self.shift}. {self.face}"

        return str_rep


class Tile:
    def __init__(self, tile_center: np.ndarray, shift_face_list: list[ShiftFace], tile_id=None) -> NoneType:
        self.shift_face_list = shift_face_list
        self.id = tile_id

        for index in range(len(self.shift_face_list)):
            self.shift_face_list[index].tiles.append(self)
    
    @staticmethod
    def build_from_shift_faces(shift_face_list: list[ShiftFace], tile_id=None) -> Tile:
        tile_center = np.mean(np.array([shift_face.get_position() for shift_face in shift_face_list]))

        new_tile = Tile(
            tile_center=tile_center,
            shift_face_list=shift_face_list,
            tile_id=tile_id,
        )

        return new_tile

    
    def __repr__(self,) -> str:
        str_rep = f"Tile {self.id}\n"
        for shift_face in self.shift_face_list:
            str_rep += str(shift_face) + '\n'
        
        str_rep = str_rep.strip()

        return str_rep


class Tiling:
    def __init__(self, name: str, basis: np.ndarray, node_list: NodeList, edge_list: EdgeList, face_list: FaceList, tiles: list[Tile]) -> NoneType:
        self.name = name
        self.basis = basis
        self.node_list = node_list
        self.edge_list = edge_list
        self.face_list = face_list
        self.tiles = tiles
    
    @staticmethod
    def build_from_json_face_dict_list(json_tiling: dict) -> Tiling:
        name = json_tiling["netName"]
        basis = Basis(
            np.array(json_tiling["unitCellBasis"]),
        )
        face_dict_list = json_tiling["faces"]

        node_list = NodeList.build_from_json(
            face_dict_list=face_dict_list,
            basis=basis,
            epsilon=1e-2
        )

        edge_list = EdgeList.build_from_json(
            face_dict_list=face_dict_list,
            basis=basis,
            node_list=node_list,
            epsilon_center=1e-2,
            epsilon_vector=1e-2,
        )

        face_list = FaceList.build_from_json(
            face_dict_list=face_dict_list,
            basis=basis,
            node_list=node_list,
            edge_list=edge_list,
            epsilon_center=1e-2,
            epsilon_vector=1e-2,
        )

        tiles_dict = {}
        for face_dict in face_dict_list:
            positions = basis.global_to_rel(np.array(face_dict["nodePositions"]))
            shift_face = ShiftFace.build_from_positions(
                positions=positions,
                node_list=node_list,
                edge_list=edge_list,
                face_list=face_list,
                basis=basis,
                epsilon_center=1e-2,
                epsilon_vector=1e-2,
            )

            tile_index = face_dict["tileIndex"]
            if tile_index not in tiles_dict:
                tiles_dict[tile_index] = []
            
            tiles_dict[tile_index].append(shift_face)
        
        tiles = [Tile.build_from_shift_faces(shift_face_list=val, tile_id=index) for (index, val) in enumerate(tiles_dict.values())]

        new_tiling = Tiling(
            name=name,
            basis=basis.vectors,
            node_list=node_list,
            edge_list=edge_list,
            face_list=face_list,
            tiles=tiles,
        )

        return new_tiling
    
    def __repr__(self,) -> str:
        str_rep = "Tiling\n"
        for tile in self.tiles:
            str_rep += str(tile) + '\n'

        str_rep = str_rep.strip()

        return str_rep



def test():
    with open('/Users/joshgoldman/Documents/Research/AI4ChemS/LQGReader/LQGReader/sample_data/aco.json', 'r') as f:
        json_tiling = json.load(f)

    tiling = Tiling.build_from_json_face_dict_list(
        json_tiling=json_tiling
    )

    print(tiling)


if __name__ == "__main__":
    test()
