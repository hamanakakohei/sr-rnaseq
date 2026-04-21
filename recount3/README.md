# recount3

ドキュメンテーション：
https://rna.recount.bio/docs/index.html

recount2のDEGなど一連の解析：
https://f1000research.com/articles/6-1558/v1
https://bioconductor.org/packages/3.21/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.R
https://bioconductor.org/packages/3.21/workflows/vignettes/recountWorkflow/inst/doc/recount-workflow.html
（2と3では、アライナーとリードカウントの方法が微妙に違うらしい（https://bioconductor.org/packages/release/bioc/manuals/recount3/man/recount3.pdf））

recount3のクイックスタートガイド：
https://bioconductor.org/packages/release/bioc/vignettes/recount3/inst/doc/recount3-quickstart.html

プロジェクトの選び方やDEG解析の準備：
https://kazumaxneo.hatenablog.com/entry/2022/02/08/024623

get_and_plot_bw.R
- bigwigをダウンロードして可視化する
- getTotalMappedで補正している


human_projectsの列の説明：
- project："<DRP or SRP or ERP>xxx"ならsraのデータ？gtexやtcgaの場合ここに30種くらいの組織名、この列でデータが一意に決まる
- organism："human"のみ
- file_source；"<gtex or sra or tcga>"
- project_home："data_sources/<gtex or sra or tcga>"
- project_type："data_sources"のみ
- n_samples：サンプル数のみ
