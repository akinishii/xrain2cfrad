# xrain2cfrad
A python script to convert XRAIN RAW/Intermidiate format to CF-Radial 1.5<br>
XRAINのRAWデータと1次処理データをCf-Radialに変換するPythonスクリプト

# 必要なもの (Requirement)
Python > 3.6  
numpy  
netCDF4  

# 使い方 (Useage)
1. P008ファイルのアーカイブ(DIASから取得できる，RAWデータのtarfile)とR005ファイルのアーカイブを(同様に取得できる，一次処理データのtarfile)を同じディレクトリに保存する。 
  
2. 以下のように実行(path/to/P008はP008アーカイブのパスに置き換えること)する
<code>python3 Conv_xrain2cfrad.py path/to/P008</code>  
  
3. "Convert success: Saved to cfrad.~.nc"(~はP008ファイルから取得した場所，日付，仰角番号情報)と表示されたら変換成功  
スクリプトを実行したディレクトリに変換されたファイルが保存されているはず  
  
# その他 (Ohters)
正常に変換されているか確認するためのコード"check_xrain_cfrad.py"を用意しました。あまりきれいな図は出ませんが確認用にどうぞ。
  
# Todo
・P008，R005のみの変換に対応させる<br>
・コードの可読性向上

# 参考
jmardr_cfradial: https://github.com/wm-ytakano/jmardr_cfradial
Cf-Radial ver1.5 document: https://github.com/NCAR/CfRadial/blob/master/docs/CfRadialDoc.v1.5.20211201.pdf
水防災オープンデータ提供サービスMPレーダ雨量計RAW・１次処理データ共通フォーマット: http://www.river.or.jp/koeki/opendata/data/03_RAWdata_Data%20format%20specifications_1.4.pdf
