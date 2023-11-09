# xrain2cfrad
A python script to convert XRAIN RAW/Intermidiate format to CF-Radial 1.5<br>
XRAINのRAWデータと1次処理データをCf-Radialに変換するPythonスクリプト

# 必要なもの
Python > 3.6  
numpy  
netCDF4  

# 使い方
1. P008ファイル(RAWデータがアーカイブされたtarfile)とR005ファイル(一次処理データ)を同じディレクトリに保存  
  
2. 以下のように実行(path/to/P008はP008ファイルのパスに置き換えること)  
<code>python3 Conv_xrain2cfrad.py path/to/P008</code>  
  
3. "Convert success: Saved to cfrad.~.nc"(~はP008ファイルから取得した場所，日付，仰角番号情報)と表示されたら変換成功  
スクリプトを実行したディレクトリに変換されたファイルが保存されているはず  
  
# その他
<li>正常に変換されているか確認するためのコード"check_xrain_cfrad.py"を用意しました。あまりきれいな図は出ませんが確認用にどうぞ  </li>
* 使用にはmatplotlibと<a href="https://arm-doe.github.io/pyart/index.html">arm-Pyart</a>が必要です  
  
