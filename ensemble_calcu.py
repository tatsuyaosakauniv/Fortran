import pandas as pd

# データの読み込み
data = pd.read_csv("flow_check_top_1104.dat", delim_whitespace=True, header=None)

# 初期設定
output = []

# 各開始行（1行目、6行目、11行目、...、86行目）から18回の処理を行う
for start_offset in range(0, 90, 5):
    current_row = start_offset  # 開始行を設定
    
    while current_row + 100 <= len(data):  # 50000行まで処理
        # tmpの設定：開始行の2列目のデータ
        tmp = data.iloc[current_row, 1]

        # tmp * 現在の行の値を1~100の番号で計算し、出力リストに追加
        for i in range(100):
            output.append([(i + 1)*0.1, tmp * data.iloc[current_row + i, 1]])

        # 次のブロックに進む（100行後に進む）
        current_row += 100

# 結果を新しいDataFrameに格納し、ファイルに書き出し
output_df = pd.DataFrame(output, columns=["Index", "Value"])
output_df.to_csv("newfile.dat", index=False, header=False)

import pandas as pd
import matplotlib.pyplot as plt

# 前のステップで作成したファイルを読み込み
data = pd.read_csv("newfile.dat", header=None, names=["Index", "Value"])

# 1列目ごとに2列目の値を合計
summed_data = data.groupby("Index")["Value"].sum().reset_index()

# 結果の確認
print(summed_data)

# 必要であれば新しいファイルとして保存
summed_data.to_csv("summed_data.csv", index=False, header=False)

# グラフの描画
plt.figure(figsize=(6, 6))
plt.plot(summed_data["Index"], summed_data["Value"], marker='', linestyle='-')
plt.xlabel("Time  ps")
plt.ylabel("J  (t) * J  (0)")

# y=0の点線を表示
plt.axhline(0, color='gray', linestyle='--')  # 点線を引く

# グリッドを消す
plt.grid(False)

# 画像として保存
plt.savefig("ensemble2.png")  # ここで画像を保存します
plt.show()