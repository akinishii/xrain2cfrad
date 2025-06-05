# xrain2cfrad
A python script to convert XRAIN RAW/Intermidiate files to CF-Radial 1.5<br>
XRAINのRAWデータと1次処理データをCf-Radialに変換するPythonスクリプト

## 必要なもの (Requirements)
- Python >= 3.6
- numpy >= 2.0
- netCDF4

---


## 使い方 (Useage)
### 単一のファイルを処理する場合
1. P008ファイルのアーカイブ(DIASから取得できる，RAWデータのtarfile)とR005ファイルのアーカイブを(同様に取得できる，一次処理データのtarfile)を同じディレクトリに保存する。 
  
2. 以下のように実行(`path/to/P008`はP008のtarのパス)する(注：R005も同じディレクトリに保存する必要あり)

    <code>python3 Conv_xrain2cfrad.py path/to/P008</code>
    
    `-o` または `--outdir` オプションで出力ファイルの保存場所を指定できる。以下の例では`./converted_nc`に保存。

    ```bash
    python3 Conv_xrain2cfrad.py path/to/P008 -o ./converted_nc
    ```
  
3. "Convert success: Saved to cfrad.~.nc"(~はP008ファイルから取得した場所，日付，仰角番号情報)と表示されたら変換成功。

    デフォルトでは`./out_nc`ディレクトリにCf-Radial変換されたファイルが保存される。

    出力ファイルのフォーマット(XRAIN田村レーダーの場合)

    `cfrad.TAMURA0000-yyyymmdd-hhmm-ELxxxxxx-DEGeee.nc`
    
    `xxxxxx`: 仰角番号, `eee`: 仰角の10倍値

### 複数のファイルを処理する場合
1. `Conv_xrain2cfrad.py`から <b>P008ファイル(.tar)</b> の相対パスを記載したテキストファイル(`filelist.txt`)を作成する

2. 以下のように実行する ※txtファイルを指定すると自動的に複数ファイル処理モードになる。

    ```bash
    python3 Conv_xrain2cfrad.py filelist.txt
    ```


## 注意 (Notes)
レーダーによっては仰角番号が同じで場合でも仰角の値が少し(0.1~0.2度)ずれる場合があります。これは、XRAINのヘッダーに収録されている仰角情報がデータによって異なることにより生じるものです。


## その他 (Ohters)
* 正常に変換されているか確認するためのコード`check_xrain_cfrad.py`を用意しました。あまりきれいな図は出ませんが確認用にどうぞ。
* REF, VEL, ZDR, KDP，RHOHV, WIDTHをそれなりにきれいに描画する`Draw_XRAIN_cfrad.py`も用意しました。

    ※上記スクリプトの実行には`arm_pyart` (https://github.com/ARM-DOE/pyart) が必要です。


## Todo

- ヘッダー情報出力モードの実装
- コードの可読性向上
- DIASから取ってきたデータをそのままのディレクトリ構造で処理できるモードの追加
- P008，R005のみの変換に対応させる


## 参考
- jmardr_cfradial: https://github.com/wm-ytakano/jmardr_cfradial
- Cf-Radial ver1.5 document: https://github.com/NCAR/CfRadial/blob/master/docs/CfRadialDoc.v1.5.20211201.pdf
- 水防災オープンデータ提供サービスMPレーダ雨量計RAW・１次処理データ共通フォーマット:
  http://www.river.or.jp/koeki/opendata/data/03_RAWdata_Data%20format%20specifications_1.4.pdf
