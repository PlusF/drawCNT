import bpy
import numpy as np
from ase.build import nanotube
from ase.geometry.analysis import Analysis


def delete_objs():
    for item in bpy.data.objects:
        bpy.data.objects.remove(item)
    for item in bpy.data.meshes:
        bpy.data.meshes.remove(item)
    for item in bpy.data.materials:
        bpy.data.materials.remove(item)
        

def create_atom(loc, radius, material):
    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16, radius=radius, calc_uvs=True, enter_editmode=False, align='WORLD', location=loc, rotation=(0.0, 0.0, 0.0), scale=(1.0, 1.0, 1.0))
    configurate_object(material)


def create_bond(loc, theta, radius, material):
    bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius=radius, depth=1.46, end_fill_type='NGON', calc_uvs=True, enter_editmode=False, align='WORLD', location=loc, rotation=theta, scale=(1.0, 1.0, 1.0))
    configurate_object(material)
    

def configurate_object(material):
    mesh = bpy.context.object.data
    for f in mesh.polygons:
        f.use_smooth = True
    mesh.materials.append(material)


def get_materials():
    carbon = bpy.data.materials.new('color')
    carbon.diffuse_color = (0.05, 0.05, 0.05, 1)
    carbon.roughness = 1.0
    boron = bpy.data.materials.new('color')
    boron.diffuse_color = (255/255, 0, 240/255, 1)
    boron.roughness = 1.0
    nitrogen = bpy.data.materials.new('color')
    nitrogen.diffuse_color = (0, 10/255, 255/255, 1)
    nitrogen.roughness = 1.0
    bond = bpy.data.materials.new('color')
    bond.diffuse_color = (0.05, 0.05, 0.05, 1)
    bond.roughness = 1.0
    
    return {
        'Carbon': carbon,
        'Boron': boron,
        'Nitrogen': nitrogen,
        'Bond': bond,
    }


def calc_rotation(coord_1, coord_2):
    if coord_1[2] < coord_2[2]:
        coord_1, coord_2 = coord_2, coord_1

    center = (coord_1 + coord_2) / 2
    length = np.linalg.norm(coord_1 - coord_2)
    theta_y = np.arccos((coord_1[2] - center[2]) / (length / 2))

    half = np.array([0, 0, length / 2])
    before = center + half

    angle = np.array([0, theta_y, 0])
    after = rotate_3d((before - center), angle) + center

    theta_loc = np.arctan2(coord_1[1] - center[1], coord_1[0] - center[0])
    theta_after = np.arctan2(after[1] - center[1], after[0] - center[0])
    theta_z = theta_loc - theta_after

    return np.array([0, theta_y, theta_z])


def rotate_3d(coord: np.array, angle: np.array):
    px, py, pz = angle

    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(px), np.sin(px)],
        [0, -np.sin(px), np.cos(px)]
    ])
    Ry = np.array([
        [np.cos(py), 0, -np.sin(py)],
        [0, 1, 0],
        [np.sin(py), 0, np.cos(py)],
    ])
    Rz = np.array([
        [np.cos(pz), np.sin(pz), 0],
        [-np.sin(pz), np.cos(pz), 0],
        [0, 0, 1],
    ])

    return Rz * Ry * Rx @ coord


def create_cnt(n, m, l):
    cnt = nanotube(n, m, length=l)
    cnt.pbc = [False, False, False]
    
    # to retrieve bond information
    ana = Analysis(cnt)
    CCBonds = ana.get_bonds('C', 'C', unique=True)[0]
    center_array = [(cnt.positions[pair[0]] + cnt.positions[pair[1]]) / 2 for pair in CCBonds]
    theta_array = [calc_rotation(cnt.positions[pair[0]], cnt.positions[pair[1]]) for pair in CCBonds]
    bond_array = np.concatenate([center_array, theta_array], axis=1)
    
    # draw atoms
    for i, loc in enumerate(cnt.positions):
        create_atom(loc, 0.17, materials['Carbon'])
    # draw bonds
    for i, (loc, theta) in enumerate(zip(center_array, theta_array)):
        create_bond(loc, theta, 0.06, materials['Bond'])


def main():
    create_cnt(6, 5, 1)
    create_cnt(7, 6, 1)
    create_cnt(8, 7, 1)
    create_cnt(9, 8, 1)


if __name__ == '__main__':
    # clear everything
    delete_objs()
    materials = get_materials()
    main()
