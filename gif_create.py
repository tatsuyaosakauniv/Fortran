from PIL import Image
import os

# BMPファイルが保存されているディレクトリを指定
bmp_directory = 'path_to_bmp_files'
gif_path = 'output.gif'

# BMPファイルを読み込む
images = []
for file_name in sorted(os.listdir(bmp_directory)):
    if file_name.endswith('.bmp'):
        file_path = os.path.join(bmp_directory, file_name)
        images.append(Image.open(file_path))

# GIFファイルを保存
if images:
    images[0].save(
        gif_path,
        save_all=True,
        append_images=images[1:],
        optimize=False,
        duration=100,
        loop=0
    )

print(f'GIFファイルが {gif_path} に保存されました')
