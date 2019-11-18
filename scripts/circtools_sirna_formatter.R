#!/usr/bin/env Rscript

# suppress loading messages
suppressMessages(library(formattable))
suppressMessages(library(kableExtra))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(colortools))

# switch to red warning color if more blast hits are found
high_count_number = 0

args <- commandArgs(trailingOnly = TRUE)
experiment_name <- args[4]
#experiment_name <- "testing"

# set output to HTML
options(knitr.table.format = 'html')

#print("script Opened")

# generic HTML header with bootstrap JS libraries
#############################################################################################################
html_header="
<html>
<head>
<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">
<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js\"></script>
<script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>
<script>
$(document).ready(function(){
$('[data-toggle=\"popover\"]').popover();
});
</script>
<script>
$(document).ready(function(){
$('[data-toggle=\"tooltip\"]').tooltip();
});
</script>
<style type=\"text/css\">
/* The max width is dependant on the container (more info below) */
.popover{
max-width: 50%; /* Max Width of the popover (depending on the container!) */
}
</style>
</head>
<body>"

html_header <- paste(html_header,"<h1>circtools siRNA design results for experiment ", experiment_name , "</h1>", sep="")

#############################################################################################################

# generate a divergent color scale with 11 shades
color_palette <- rev(brewer.pal(n = 5, name = 'RdBu'))

#default TM value
default_score <- 50
default_seed_stability <- 21.5
default_thermo_stability <- 0

#construct_color_column method should give a color based on input threshholds
construct_color_column <- function(column, default_value, palette, comp){
  color_column = c()
  if(comp == 1){
    for(val in column){
      if(val > default_value){
        color_column <- append(color_column, "#ff726f")
      }
      else{
        color_column <- append(color_column, "#32CD32")
      }
    } 
  } 
  if(comp == 0){
    for(val in column){
      if(val < default_value){
        color_column <- append(color_column, "#ff726f")
      }
      else{
        color_column <- append(color_column, "#32CD32")
      }
    }
  }
  return(as.character(color_column))
}

construct_color_column2 <- function(column, default_value, palette)
{
  top_val <- (max(column, na.rm=T) - default_value)
  bottom_val <- (default_value - min(column, na.rm=T))
  
  if (top_val > bottom_val){
    from <- default_value - top_val
    to <- default_value + top_val
  } else if (top_val == bottom_val) {
    from <- default_value - 1
    to <- default_value + 1
  } else {
    from <- default_value - bottom_val
    to <- default_value + bottom_val
  }
  return(as.character(cut(column,seq( from, to, length.out= length(palette)+1 ), labels=palette, include.lowest = T)))
}

construct_color_column3 <- function(column)
{
  color_column = c()
  for(val in column){
    color_column <- append(color_column, "#ffffff")
  }
  return(as.character(color_column))
}

generate_complementary_column <- function (input){
  return (complementary(input)[2])
}

# read data file name from args

#args <- commandArgs(trailingOnly = TRUE)
csv <- args[1]
blast_r <- args[2]
output_dir <- args[3]
experiment_title <- args[4]
#csv <- file("stdin", 'r')
#open(args[2])

# read whole file into data table
data_table <- read.csv(csv, skip=1, header = FALSE)#, sep = "\t")


data_table$circid <- paste(data_table$V2,data_table$V3,data_table$V4,data_table$V5,data_table$V6,sep="_")
data_table$circid <- paste(sep="","<img src=",data_table$circid,data_table$V8,".svg>")

#write(as.character(data_table$V7), file = "viewing_file")
#print(data_table)
# remove unused columns
write(as.character(data_table$circid), file = "viewing_file")
data_table <- data_table[-c(1,8)]



# correctly name the columns
colnames(data_table) <- c(  "Annotation",
                            "Chr",
                            "Start",
                            "Stop",
                            "Strand",
                            "siRNA",
                            "Silencing_Score",
                            "Rule",
                            "Blast",
                            "Seed_Duplex_Stability",
                            "Thermodynamic_Stability",
                            "ID"
)


#data_table$score  = construct_color_column(as.numeric(as.character(data_table$Silencing_Score)),default_score,color_palette, 0)
#data_table$left_tm_color   = construct_color_column(data_table$TM_left,default_tm_value,color_palette)

#data_table$seedStability   = construct_color_column(as.numeric(as.character(data_table$Seed_Duplex_Stability)),default_seed_stability,color_palette, 1)
#data_table$right_gc_color  = construct_color_column(data_table$GC_right,default_gc_value,color_palette)

#data_table$thermoStability  = construct_color_column(as.numeric(data_table$Thermodynamic_Stability),default_thermo_stability,color_palette, 1)

if(identical(blast_r, "True")){
    data_table$blast  = construct_color_column2(as.numeric(as.character(data_table$Blast)),default_thermo_stability,color_palette)#, 1)
}else{
    data_table$blast = construct_color_column3(as.character(data_table$Blast))
}
 
colnames_final <- c(        "Annotation",
                            "Chr",
                            "Start",
                            "Stop",
                            "Strand",
                            "siRNA",
                            "Silencing Score",
                            "Rule",
                            "Blast",
                            "Seed-Duplex Stability",
                            "Thermodynamic Stability"
                            #"Forward",
                            #"BLAST",
                            #"Reverse",
                            #"BLAST"
)

# run_primer_design a column with BLAST hit counts
#data_table$BLAST_left_count <- lengths(regmatches(data_table$BLAST_left, gregexpr(";", data_table$BLAST_left))) + 1
#data_table$BLAST_right_count <- lengths(regmatches(data_table$BLAST_right, gregexpr(";", data_table$BLAST_right))) + 1

#data_table$BLAST_left_count[data_table$BLAST_left_count == 1] = 0
#data_table$BLAST_right_count[data_table$BLAST_right_count == 1] = 0

# replace ; with HTML linebreaks for hover popover text
#data_table$BLAST_left <- gsub(";", "<br/><br/>", data_table$BLAST_left)
#data_table$BLAST_right <- gsub(";", "<br/><br/>", data_table$BLAST_right)

# remove 0 entries from location columns for provided circRNA FASTA files
data_table$Chr <- gsub("\\b0\\b", "", data_table$Chr )
data_table$Start <- gsub("\\b0\\b", "", data_table$Start )

data_table$Stop  <- ifelse( data_table$Start == "" , "", data_table$Stop )

data_table$Strand<- gsub("\\b0\\b", "", data_table$Strand )

# clean up rownames to hide them lateron
rownames(data_table) <- c()


# main output table generation
output_table <- data_table %>%
  
  mutate(
    Silencing_Score = cell_spec(Silencing_Score, background= ifelse(Silencing_Score <= 49, "#ff726f", "#32CD32")),
    siRNA = cell_spec(siRNA, escape = F, popover = spec_popover( title = "Graphical represensation of siRNA\"data-html=\"True\"", position = "right", content =ID )),
    #Forward = cell_spec(escape = F, Left_, popover = spec_popover( title = "Graphical represensation of designed primers and annotated circRNA structure\"data-html=\"True\"", position = "left", content =ID ), background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
    #                    color = ifelse(BLAST_left_count > high_count_number, "white", "white")),
    
    #L = cell_spec(paste(BLAST_left_count),
    #              popover = spec_popover(content = BLAST_left, title = "Blast Hits\"data-html=\"True\"", position = "right"),
    #              background = ifelse(BLAST_left_count > high_count_number, "red", "darkgreen"),
    #              color = ifelse(BLAST_left_count > high_count_number, "white", "white"), bold = "true"),
    
    #Reverse = cell_spec(escape = F, Right_, popover = spec_popover( title = "Graphical represensation of designed primers and annotated circRNA structure\"data-html=\"True\"", position = "left", content =ID ), background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
    #                    color = ifelse(BLAST_right_count > high_count_number, "white", "white")),
    
    #R = cell_spec(paste(BLAST_right_count),
    #              popover = spec_popover(content = BLAST_right, title = "Blast Hits\"data-html=\"True\"", position = "left"),
    #              background = ifelse(BLAST_right_count > high_count_number, "red", "darkgreen"),
    #              color = ifelse(BLAST_right_count > high_count_number, "white", "white"), bold = "true"),
    
    Seed_Duplex_Stability = cell_spec(Seed_Duplex_Stability, background= ifelse(Seed_Duplex_Stability >= 21.5, "#ff726f", "#32CD32")),
    #Thermodynamic_Stability = color_bar(thermoStability)(Thermodynamic_Stability),
    Thermodynamic_Stability = cell_spec(Thermodynamic_Stability, background= ifelse(Thermodynamic_Stability > 0, "#ff726f", "#32CD32")),
    #Thermodynamic_Stability = fixNegativeSigns(Thermodynamic_Stability),
    
    
    
    Blast = color_bar(blast)(Blast),
    
    #GC_left = color_bar(left_gc_color)(GC_left),
    #GC_right = color_bar(right_gc_color)(GC_right),
    
    
    Strand = formatter('span', style = style(font.weight = "bold"))(Strand)
  ) %>% 
  #select(- siRNA) %>%
  #select(- score) %>%
  #select(- Rule) %>%
  #select(- seedStability) %>%
  #select(- thermoStability) %>%
  select(- blast) %>%
  select(- ID) %>%
  select(Annotation, everything()) %>%
  kable("html", escape = F, col.names=colnames_final, align=c(rep('l',times=4), 'c', 'r', 'r', 'l', 'r', 'r', 'r')) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = T) 

fixNegativeSigns <- function(vector){
  newVector <- c()
  for(i in seq(0,length(vector),1)){
    x = vector[i]
    if(identical(substr(x,nchar(x),nchar(x)),"-")){
      x = paste("-",substr(x,0,nchar(x)-1), sep="")
      newVector <- c(newVector, x)
    }
  }
  return(newVector)
}




#kable_styling(position='float_right',full_width=F)
# column_spec(5, width = "3cm")
#add_header_above(c("Input circRNAs" = 5, "Designed Primers" = 9)) # %>%
# group_rows("Group 1", 4, 7) %>%
# group_rows("Group 1", 8, 10)
# collapse_rows(columns = 1)

#setwd("/beegfs/")
#print(data_table)
outputString = paste(html_header, output_table, sep="")
#print(outputString)
html_file = paste(output_dir,experiment_title, ".html", sep = "")
write(paste(html_header, output_table, sep=""), file = html_file)
