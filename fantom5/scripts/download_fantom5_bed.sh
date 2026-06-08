#!/bin/bash
# FANTOM5 hg38 basic BED ファイル一括ダウンロードスクリプト
# 使い方: bash download_fantom5_bed.sh [保存先ディレクトリ]
#   例: bash download_fantom5_bed.sh ./fantom5_bed

set -euo pipefail

BASE_URL="https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/basic"

DIRS=(
    "human.cell_line.LQhCAGE"
    "human.cell_line.hCAGE"
    "human.fractionation.hCAGE"
    "human.primary_cell.LQhCAGE"
    "human.primary_cell.hCAGE"
    "human.timecourse.LQhCAGE"
    "human.timecourse.hCAGE"
    "human.tissue.hCAGE"
)

OUT_DIR="${1:-./fantom5_bed}"
URL_LIST=$(mktemp /tmp/aria2_fantom5_XXXXXX.txt)

echo "=== FANTOM5 BED ダウンロード ==="
echo "保存先: ${OUT_DIR}"
echo ""

# Python でディレクトリ一覧を解析し aria2c 用リストを生成
python3 - "${OUT_DIR}" "${URL_LIST}" "${BASE_URL}" "${DIRS[@]}" <<'PYEOF'
import sys, re, urllib.request, urllib.parse

out_dir   = sys.argv[1]
list_file = sys.argv[2]
base_url  = sys.argv[3]
dirs      = sys.argv[4:]

def sanitize(name):
    """スペース・記号類を _ に置換してファイル名を安全にする"""
    # まずURLデコード
    name = urllib.parse.unquote(name)
    # 記号をアンダースコアに置換 (英数字・ドット・ハイフン・_ 以外)
    name = re.sub(r'[^\w.\-]', '_', name)
    # 連続する _ をまとめる
    name = re.sub(r'_+', '_', name)
    return name

entries = []
for d in dirs:
    dir_url = f"{base_url}/{d}/"
    print(f"[リスト取得] {dir_url}", flush=True)
    try:
        with urllib.request.urlopen(dir_url, timeout=30) as r:
            html = r.read().decode("utf-8", errors="replace")
    except Exception as e:
        print(f"  WARNING: 取得失敗 ({e})", flush=True)
        continue

    # href属性からファイル名を抽出
    hrefs = re.findall(r'href="([^"]+\.bed\.gz)"', html)
    print(f"  -> {len(hrefs)} ファイル", flush=True)
    for href in hrefs:
        # href の値(2重エンコード済み)をそのままURLに使う
        # Apache が受け取って1回デコードし、実ファイル名(1重エンコード)にマッピングする
        full_url = dir_url + href
        # ローカルファイル名: 2回デコードして平文にしてからサニタイズ
        plain = urllib.parse.unquote(urllib.parse.unquote(href))
        local_name = sanitize(plain)
        entries.append((full_url, f"{out_dir}/{d}", local_name))

with open(list_file, "w") as f:
    for url, save_dir, fname in entries:
        f.write(f"{url}\n")
        f.write(f"  dir={save_dir}\n")
        f.write(f"  out={fname}\n")

print(f"\n合計: {len(entries)} ファイル")
PYEOF

echo ""
echo "[ダウンロード開始]"
aria2c \
    --input-file="${URL_LIST}" \
    --max-concurrent-downloads=16 \
    --split=1 \
    --max-connection-per-server=2 \
    --retry-wait=5 \
    --max-tries=5 \
    --auto-file-renaming=false \
    --lowest-speed-limit=50K \
    --console-log-level=notice \
    --summary-interval=60

rm -f "${URL_LIST}"
echo ""
echo "=== 完了 ==="
