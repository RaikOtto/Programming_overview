library("ggplot2")
library("RMySQL")
library("googleVis")

stringsAsFactors=F

### get data

con=dbConnect(MySQL(),user="jevgeni",password="test",dbname="gender_vis",host="textmining")

# main data
rs_main_data_raw =dbSendQuery(con, "SELECT * FROM gender_vis.author_authorship_comp;")
main_data_raw=fetch(rs_main_data_raw,n=10000,header=T)

# timeline data
rs_motion_data =dbSendQuery(con, "SELECT * FROM gender_vis.frac_contrib_country_remap_non_china_per_pub_year;")
motion_data_raw=data.frame(fetch(rs_motion_data,n=10000,header=T),stringsAsFactors=F)
motion_data_raw$pub_year = as.numeric(motion_data_raw$pub_year)

# difference data

rs_diff_data =dbSendQuery(con, "SELECT * FROM gender_vis.authornet_nature_comp;")
diff_data_raw=data.frame(fetch(rs_diff_data,n=10000,header=T),stringsAsFactors=F)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
diff_data_raw= data.frame(apply(diff_data_raw,2,FUN=trim),stringsAsFactors=F)
diff_data_raw$authornet_fmratio = as.numeric( diff_data_raw$authornet_fmratio )
diff_data_raw$nature_fmratio = as.numeric( diff_data_raw$nature_fmratio )
diff_data_raw$difference = as.numeric( diff_data_raw$difference )

# count data

rs_count_data  = dbSendQuery(con, "SELECT * FROM gender_vis.country_gender_frac_country_remap;")
count_data_raw = data.frame(fetch(rs_count_data,n=10000,header=T),stringsAsFactors=F)
count_data_raw = count_data_raw[order(count_data_raw$sum_count,decreasing=T),]

# close connection 
dbDisconnect(con)

### preprocess data

colnames(main_data_raw) = c("Country","FMR_Author","FMR_Authorship","Amount_Authors","Amount_Authorships","Google_Code")

Dif_Numb_FMR    = round( as.numeric(main_data_raw$FMR_Author) - as.numeric(main_data_raw$FMR_Authorship),2)

main_data_raw = data.frame( stringsAsFactors=F,
                as.character(main_data_raw$Country),
                as.character(main_data_raw$Google_Code),
                main_data_raw$FMR_Author,
                main_data_raw$FMR_Authorship,
                Dif_Numb_FMR,
                main_data_raw$Amount_Authors,
                main_data_raw$Amount_Authorships
)
        
colnames(main_data_raw) = c(
        "Country",
        "Country_Code",
        "FMR_Author",
        "FMR_Authorship",
        "FMR_Difference",
        "Amount_Authors",
        "Amount_Authorships"
)

main_data_raw$FMR_Author                        = as.numeric(main_data_raw$FMR_Author)
main_data_raw$FMR_Authorship            = as.numeric(main_data_raw$FMR_Authorship)
main_data_raw$Dif_Numb_FMR                      = round( as.numeric(main_data_raw$FMR_Author) - as.numeric(main_data_raw$FMR_Authorship),2)
main_data_raw$Amount_Authors            = as.numeric(main_data_raw$Amount_Authors)
main_data_raw$Amount_Authorships        = as.numeric(main_data_raw$Amount_Authorships)

### start the stuff

shinyServer(function(input, output, session) {

        observe({
                
                # just filtering 
                main_data       = main_data_raw#[ main_data_raw$Amount_Authors >= input$num_pubs,]
                diff_data       = diff_data_raw[ diff_data_raw$country %in% main_data$Country,]
                #write.table(diff_data_raw,"test.tab")
                motion_data = motion_data_raw#[ as.integer(motion_data_raw$num_authors )>=input$num_pubs,]
                count_data      = count_data_raw[ count_data_raw$country %in% main_data$Country, ]

                main_data       = main_data[ order(as.numeric(main_data$Amount_Authors),decreasing=T),]
                
                # type casts, to be changed later
                
                main_data_raw$FMR_Author                = as.numeric(main_data_raw$FMR_Author)
                main_data_raw$FMR_Authorship    = as.numeric(main_data_raw$FMR_Authorship)
                main_data_raw$Dif_Numb_FMR              = round( as.numeric(main_data_raw$FMR_Author) - as.numeric(main_data_raw$FMR_Authorship),2)
                main_data_raw$Amount_Authors    = as.numeric(main_data_raw$Amount_Authors)
                main_data_raw$Amount_Authorships= as.numeric(main_data_raw$Amount_Authorships)
                
                # here we go
                
                #if (input$scale_bar == "Log10"){
                
                 #      main_data$Amount_Authors                        = round(log10( as.integer(main_data$Amount_Authors)),2)
                 #     main_data$Amount_Authorships            = round(log10( as.integer(main_data$Amount_Authorships)),2)
                #}
                
                #vis_data = switch(input$radio_buttons,

                 #              "Counts F/M Researchers per country our method"= with(main_data,data.frame(cbind( 
                 #                      Country,Amount_Authors, Amount_Authorships))),
                 #              "FMR per country our method" = with(main_data,data.frame(cbind( 
                 #                      Country,FMR_Author,FMR_Authorship,Dif_Numb_FMR))),
                 #              "Data difference to other pub"=with(diff_data,data.frame(cbind(
                 #                      countries,main_data)))
                #)
                
                output$barPlot_Authors = renderGvis({                   
                        
                        author_data = data.frame(main_data$Country,main_data$Amount_Authors,main_data$Amount_Authorships)
                        colnames(author_data) = c("Country","Authors","Authorships")
                        gvisBarChart( author_data, xvar="Country")
                })
                
                output$barPlot_FMR = renderGvis({

                        author_data = data.frame(main_data$Country,main_data$FMR_Author,main_data$FMR_Authorship)
                        colnames(author_data) = c("Country","FMR_Author","FMR_Authorship")
                        gvisBarChart( author_data, xvar="Country")
                })

                output$barPlot_gender = renderGvis({

                        nr_data = data.frame(count_data$country,count_data$male_count,count_data$female_count)
                        colnames(nr_data) = c("Country","Male","Female")
                        gvisBarChart( nr_data, xvar="Country")
                })

                output$barPlot_Differences = renderGvis({

                        vis_data = diff_data[,seq(1,3)]
                        #vis_data$difference = abs(vis_data$difference)
                        vis_data = vis_data[order(vis_data$authornet_fmratio,decreasing=T),]

                        gvisBarChart( vis_data, xvar="country",
                        options = list( width=800, height=800 ) )
                })
                
                output$motionPlot = renderGvis({

                        gvisMotionChart( motion_data, timevar="pub_year", idvar="country",
                        options = list( width=600, height=400 )
                        )
                })
                
                output$scatter_plot = renderPlot({
                
                        qp = qplot( data=main_data, x=FMR_Authorship,y= FMR_Author,size=Amount_Authors)
                        qp = qp + geom_smooth() + geom_abline(col="red", intercept = .35)
                        qp = qp + ggtitle("Correlation between amount authors and authorship for selected countries")
                        print(qp)
                })
                
                output$world_map = renderGvis({
                        
                        gvisGeoChart(
                                data=main_data,
                                colorvar="FMR_Author",
                                locationvar="Country_Code",
                                options=list(
                                #title="Hello World",
                                #titleTextStyle="{color:red,fontName:Courier,fontSize:16}",
                                colorAxis="{values:[.5,1.5],colors:['blue','red']}",
                                #title="FMR per country",
                                width=800, height=800,align="top")
                                )
                })
                
                output$world_map_authorships = renderGvis({
                        
                        gvisGeoChart(
                                data=main_data,
                                colorvar="FMR_Authorship",
                                locationvar="Country_Code",
                                options=list(colorAxis="{values:[.15,.94],colors:['blue','red']}",
                                title="FMR Authorhips per country",
                                width=800, height=800,align="top")
                        )
                })
                
                output$world_map_difference_nature_us = renderGvis({
                        
                        gvisGeoChart(
                                data=diff_data,
                                colorvar="difference",
                                locationvar="code",
                                options=list(
                                width=800, height=800)
                                )
                })
                
                output$table_data = renderDataTable({
                        vis_data
                })
                
                output$count_table_data = renderDataTable({
                        vis_data = count_data[,seq(1,6)]
                })

                output$data_dif_nature_us = renderDataTable({
                        
                        vis_data = diff_data[,seq(1,4)]
                        vis_data = vis_data[order(vis_data$difference,decreasing=T),]
                        vis_data
                })
                

                # download stuff
                output$main_data_raw = downloadHandler( 
                        
                        filename = function(){"data.tab"},
                        content = function(file)( write.table( main_data_raw,file,sep="\t",row.names=F) )
                )
                
                output$time_series_data = downloadHandler( 
                        
                        filename = function(){"time_series_data.tab"},
                        content  = function(file)( write.table( motion_data,file,sep="\t",row.names=F) )
                )
                
                output$down_dif_data = downloadHandler( 
                        
                        filename = function(){"diff_data.tab"},
                        content  = function(file)( write.table( diff_data,file,sep="\t",row.names=F) )
                )
                
        })
})
