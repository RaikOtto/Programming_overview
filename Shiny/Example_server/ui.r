library("shiny")

shinyUI(
	navbarPage("Global gender (in)equality among scientists",
 
		tabPanel("Home",
			div("Here, we present a study of the gender analysis of scientists working in the 
				biomedical field. The analysis is based on publications available through PubMed and the 
				consideration of authorships as well as individual authors in countries around the world. 
				Three different aspects of the gender distribution among scientists are investigated – 
				the proportion of female and male authorships and authors in each country, the evolution 
				of the gender proportions over time as well as the number of citations male and female 
				authors receive. Details about the data and the methods are described in the FAQ. If you 
				have further questions, please do not hesitate to contact us.",align="center"
			)
		),
		navbarMenu("Author gender analysis",

				tabPanel("Comparison of author gender ratios",
					
					div(
						"Bar plot - Comparison of author gender ratios",
						htmlOutput("barPlot_Authors"),align="center"),
						sliderInput("num_pubs", min =10^3, max = 10^5,value = 10^5,
							label="Include countries with less authors"
						)
				),
				tabPanel("Comparison of author gender ratios",
					div(
						"Geo chart - Comparison of author gender ratios",
						htmlOutput("barPlot_Authors"),align="center"),
						sliderInput("num_pubs", min =10^3, max = 10^5,value = 10^5,
							label="Include countries with less authors"
					)
				),
				tabPanel("Bar plot - Comparison of authorship and author gender ratios"),
				tabPanel("Geo chart - Author gender ratio"),
				tabPanel("Geo chart - Authorship gender ratio"),
				tabPanel("Geo chart - Comparson of authorship and author gender ratios"),
				tabPanel("Scatter plot - Correlation between author and authorship gender ratios",
				htmlOutput("scatter_plot")),
				tabPanel("Bar plot -Number of authors")				
		),
		
		navbarMenu("Time course analysis",
				tabPanel("Sub-Component A"),
				tabPanel("Sub-Component B")
		),
		
		navbarMenu("Citations analysis",
				tabPanel("Sub-Component A"),
				tabPanel("Sub-Component B")
		),
		
		tabPanel("Data Download"),
		
		tabPanel("FAQ"),
		
		tabPanel("Contact")
		
	#, footer = "test"
	)
)
