# download xross correlation values
from pathlib import Path
from one.webclient import http_download_file
url = "https://mouse.brain-map.org/agea/data/P56/download?seed=6800,4200,5600&age=P56"
zip_name = "6800_4200_5600.zip"
folder_download = Path("/datadisk/gdrive/2022/08_gene_expression_atlas")

fn = http_download_file(url, target_dir=folder_download)
Path(fn).rename(zip_name)
