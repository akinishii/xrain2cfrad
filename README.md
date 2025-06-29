# xrain2cfrad
A python script to convert XRAIN RAW/Intermidiate files to CF-Radial 1.5<br>
XRAINのRAWデータと1次処理データをCf-Radialに変換するPythonスクリプト

## 必要なもの (Requirements)
- Python >= 3.6
- numpy >= 2.0
- netCDF4

---


## 使い方 (Useage)
### レーダー地点名と時間範囲を指定して実行(推奨)
run_xrain2cfrad_orgdir.pyを実行することで、DIASから取得したXRAINディレクトリからそのまま変換可能です。
1. `run_xrain2cfrad_orgdir.py`と`Conv_xrain2cfrad.py`を、DIASからDLした`XRAIN`ディレクトリと同じディレクトリに置く。

2. 変換対象のレーダー地点(DIASのDLページ記載のアルファベット名。複数指定可)と、日時の範囲を確認する。

3. `run_xrain2cfrad_orgdir.py`を以下のように実行。`start_date`や`end_date`はyyyymmddHHMM形式のJST。

   <code> python3 run_xrain2cfrad_orgdir.py SITE [SITE ...] -d start_date end_date </code>

   例) 船橋局と新横浜局の、2025年3月3日09:00～10:00 JSTのデータを変換する場合。

   <code> python3 run_xrain2cfrad_orgdir.py FUNABASHI SHINYOKO -d 202503030900 202503031000 </code>

   `-o` または `--outdir` オプションで出力ファイルの保存場所を指定できる。また、`-i`または`--inputpath`で任意のディレクトリにあるXRAINディレクトリを読み出し可能。
   
   以下の例は`/data/orgdata`にある`XRAIN`ディレクトリを読み込み、`./converted_nc`に保存する場合。
   
   ```bash
   python3 run_xrain2cfrad_orgdir.py FUNABASHI -d 202503030900 202503031000 -i /data/orgdata -o ./converted_nc
   ```

4. "Convert success: Saved to cfrad.~.nc"と表示されたファイルは変換成功。

    デフォルトでは`./out_nc`ディレクトリにCf-Radial変換されたファイルが保存される。

    出力ファイルのフォーマット(XRAIN船橋局の場合)。<b>ファイル名の時刻はJSTである点に注意！</b>

    `cfrad.FUNABASHI0-yyyymmdd-hhmm-ELxxxxxx-DEGeee.nc`
    


### 単一のファイルを処理する場合
1. P008ファイルのアーカイブ(DIASから取得できる，RAWデータのtarfile)とR005ファイルのアーカイブを(同様に取得できる，一次処理データのtarfile)を同じディレクトリに保存する。 
  
2. 以下のように実行(`path/to/P008`はP008のtarのパス)する(注：R005も同じディレクトリに保存する必要あり)

    <code>python3 Conv_xrain2cfrad.py path/to/P008</code>
    
    `-o` または `--outdir` オプションで出力ファイルの保存場所を指定できる。以下の例では`./converted_nc`に保存。

    ```bash
    python3 Conv_xrain2cfrad.py path/to/P008 -o ./converted_nc
    ```
  

### 複数のファイルを処理する場合
1. `Conv_xrain2cfrad.py`から <b>P008ファイル(.tar)</b> の相対パスを記載したテキストファイル(`filelist.txt`)を作成する

2. 以下のように実行する ※txtファイルを指定すると自動的に複数ファイル処理モードになる。

    ```bash
    python3 Conv_xrain2cfrad.py filelist.txt
    ```


## 注意 (Notes)
レーダーによっては仰角番号が同じ場合でも仰角の値が若干(0.1~0.2度)ずれる場合があります。これは、XRAINのヘッダーに収録されている仰角代表値がデータにより異なることで生じたものです。現状、仰角の代表値が不明なレーダーもあるため仰角をそろえる処理はしていません。


## その他 (Ohters)
* 正常に変換されているか確認するためのコード`check_xrain_cfrad.py`を用意しました。あまりきれいな図は出ませんが確認用にどうぞ。
* REF, VEL, ZDR, KDP，RHOHV, WIDTHをそれなりにきれいに描画する`Draw_XRAIN_cfrad.py`も用意しました。

    ※上記スクリプトの実行には`arm_pyart` (https://github.com/ARM-DOE/pyart) が必要です。


## Todo
- ヘッダー情報出力スクリプトの実装
- 変換ファイルの出力ディレクトリを柔軟にする(例: 地点別に分ける)


## 参考
- jmardr_cfradial: https://github.com/wm-ytakano/jmardr_cfradial
- Cf-Radial ver1.5 document: https://github.com/NCAR/CfRadial/blob/master/docs/CfRadialDoc.v1.5.20211201.pdf
- 水防災オープンデータ提供サービスMPレーダ雨量計RAW・１次処理データ共通フォーマット:
  http://www.river.or.jp/koeki/opendata/data/03_RAWdata_Data%20format%20specifications_1.4.pdf
