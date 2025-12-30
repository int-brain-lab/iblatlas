# %%
import urllib.request

import numpy as np
import tqdm
import trimesh

from one import params
import iblatlas.atlas


CCF_URL = ('https://download.alleninstitute.org/informatics-archive/'
           'current-release/mouse_ccf/annotation/ccf_2017/structure_meshes/')
OBJ_PATH = params.get_cache_dir().joinpath('histology', 'ATLAS', 'obj')
OBJ_PATH.mkdir(parents=True, exist_ok=True)


def dl(atlas_id):
    mesh_url = CCF_URL + str(atlas_id) + '.obj'
    fn = OBJ_PATH.joinpath(f"{atlas_id}.obj")
    try:
        urllib.request.urlretrieve(mesh_url, fn)
    except Exception as e:
        print(f"Error: {str(e)}")


def get_color(atlas_id, alpha=255, br=None):
    _, idx = br.id2index(atlas_id)
    color = br.rgb[idx[0][0], :]
    return np.hstack((color, [alpha])).astype(np.uint8)


# %%
GLB_PATH = params.get_cache_dir().joinpath(
    'histology', 'ATLAS', 'Needles', 'Allen', 'meshes.glb')
atlas = iblatlas.atlas.AllenAtlas()
scene = trimesh.Scene()

for aid in tqdm.tqdm(np.unique(np.abs(atlas.regions.id))):
    file_obj = OBJ_PATH.joinpath(f"{aid}.obj")
    if not file_obj.exists():
        continue
    mesh = trimesh.load(file_obj)
    color = get_color(aid, br=atlas.regions)
    mesh.visual.vertex_colors = color
    scene.add_geometry(mesh, node_name=f"mesh_{aid}")

scene.export(GLB_PATH)
