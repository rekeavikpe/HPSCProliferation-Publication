module ManuscriptStyle
using CairoMakie, Ghostscript_jll
include("../src/PublicationModels.jl")
using .PublicationModels
export publication_model, storage_model, model_label, publication_order
export manuscript_theme!, export_figure, COLORS, MARKERS, MODEL_COLORS, SCALE_COLOR
const SCALE_COLOR="#005A70"
const COLORS=Dict(1=>"#0072B2",2=>"#E69F00",3=>"#009E73",4=>"#CC79A7",5=>"#56B4E9",6=>"#D55E00")
const MARKERS=Dict(1=>:circle,2=>:rect,3=>:diamond,4=>:hexagon,5=>:cross,6=>:xcross)
const MODEL_COLORS=[:steelblue,:red,:orange,:goldenrod,:purple,:green,:brown]
function manuscript_theme!()
    macfont="/System/Library/Fonts/Supplemental/Arial.ttf"
    macbold="/System/Library/Fonts/Supplemental/Arial Bold.ttf"
    bundled=normpath(joinpath(dirname(pathof(CairoMakie.Makie)),"..","assets","fonts"))
    font=get(ENV,"HPSC_FONT",isfile(macfont) ? macfont : joinpath(bundled,"TeXGyreHerosMakie-Regular.otf"))
    boldfont=get(ENV,"HPSC_BOLD_FONT",isfile(macbold) ? macbold : joinpath(bundled,"TeXGyreHerosMakie-Bold.otf"))
    set_theme!(Theme(font=font,fonts=(regular=font,bold=boldfont,italic=font,bold_italic=boldfont),fontsize=8,Axis=(xlabelsize=8,ylabelsize=8,titlesize=9,titlefont=:bold,xticklabelsize=7,yticklabelsize=7,spinewidth=.7,xgridvisible=true,ygridvisible=true,xgridcolor=(:gray,.12),ygridcolor=(:gray,.12),xgridwidth=.35,ygridwidth=.35,xticksize=3,yticksize=3),Lines=(linewidth=1.3,),Scatter=(markersize=4,),Legend=(framevisible=true,framewidth=.6,labelsize=7,padding=(6,6,4,4))))
end
function export_figure(out,stem,fig;dpi=600)
    mkpath(out)
    mktempdir() do staging
        pdf=joinpath(staging,stem*".pdf")
        save(pdf,fig;pt_per_unit=1)
        save(joinpath(staging,stem*".svg"),fig;pt_per_unit=1)
        gs=Ghostscript_jll.gs()
        for (extension,device) in (("png","png16m"),("jpeg","jpeg"),("tiff","tiff24nc"),("eps","eps2write"))
            dest=joinpath(staging,stem*"."*extension)
            run(`$gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=$device -r$dpi -dJPEGQ=98 -sOutputFile=$dest $pdf`)
        end
        for file in readdir(staging;join=true)
            mv(file,joinpath(out,basename(file));force=true)
        end
    end
end
end
