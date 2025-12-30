# uv pip install trimesh
# uv pip install datoviz
import numpy as np
import trimesh

import datoviz as dvz
from iblatlas.atlas import AllenAtlas
import one.remote.aws as aws

ba = AllenAtlas()


def download_glb_file():
    file_path = AllenAtlas._get_cache_dir().joinpath('meshes.glb')
    if not file_path.exists():
        file_path.parent.mkdir(exist_ok=True, parents=True)
        aws.s3_download_file(f'atlas/{file_path.name}', file_path)
    return file_path


file_glb = download_glb_file()
scene = trimesh.load_scene(file_glb)


# %%
# Setup the app
app = dvz.App(background='white')
figure = app.figure(gui=True)
panel = figure.panel()
arcball = panel.arcball()

# Load mesh for VIS, HPF, TH and root
region_idx = [669, 1089, 549, 997]
for aid in region_idx:
    mesh = scene.geometry[f'{aid}.obj']
    mesh_pos = np.ascontiguousarray(mesh.vertices, dtype=np.float32)
    mesh_idx = np.ascontiguousarray(mesh.faces.ravel(), dtype=np.uint32)
    mesh_color = np.ascontiguousarray(mesh.visual.vertex_colors, dtype=np.uint8)
    # Convert mesh coordinates to IBL coordinate space
    mesh_pos = ba.ccf2xyz(mesh_pos, ccf_order='apdvml')
    mesh_pos = np.ascontiguousarray(mesh_pos, dtype=np.float32)
    # Add transparency to the mesh
    mesh_idx = mesh_idx.astype(np.uint32).ravel()
    mesh_color[:, 3] = 64

    # Offset both mesh and clusters so they display nicely on the window
    if aid == 997:
        mesh_color[:, 3] = 32
    mesh_pos *= 200
    # Add mesh to the display
    visual = app.mesh(indexed=True, lighting=True, cull='back')
    visual.set_data(
        position=mesh_pos,
        color=mesh_color,
        index=mesh_idx,
        compute_normals=True,
    )
    panel.add(visual)

app.run()
app.destroy()
