<?php
    $title = "CRISPOR - TEFOR&#39;s CRISPR ONLINE TOOL";
    $metadescription = "Design your Cas9 RNA guide and control the off-targets. http://tefor.net/crispor
        CRISPOR - CRISPr selectOR - is a program that helps design and evaluate target sites for use with the CRISPR/Cas9 system.
        It uses the BWA algorithm to identify guide RNA sequences for CRISPR mediated genome editing. It searches for off-target sites (with and without mismatches), shows them in a table and annotates them with flanking genes.
        CRISPOR currently includes more than 70 genomes ! To add your genome of interest, contact the CRISPOR web site manager: penigault@tefor.net.";        
        
    echo "<title>$title</title>";
    echo "<meta name='description' content='$metadescription'/>
            ";
    echo "<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />
            ";
    echo "<meta property='fb:admins' content='692090743' />
            ";
    echo "<meta property='og:title' content='$title' />
            ";
    echo "<meta property='og:description' content='$metadescription'/>
            ";
    echo "<meta property='og:type' content='website' />
            ";
    echo "<meta property='og:url' content='http://tefor.net/crispor/crispor.cgi' />
            ";
    echo "<meta property='og:image' content='http://tefor.net/crispor/image/CRISPOR.png' />
            ";
    
    echo "<script src='https://ajax.googleapis.com/ajax/libs/jquery/1.11.0/jquery.min.js'></script>
         ";
    echo "<script src='http://code.jquery.com/ui/1.11.1/jquery-ui.min.js'></script>
         ";   
    echo "<script src='https://apis.google.com/js/platform.js' async defer></script>
         ";

    
    $include[] = "http://tefor.net/main/style/style_general.css";    
    $include[] = "http://tefor.net/main/specific/style/squeleton.css";    
    $include[] = "http://tefor.net/main/style/menu.css";
    $include[] = "http://tefor.net/main/specific/style/menu.css";
    $include[] = "http://tefor.net/main/specific/style/logo.css";
//    $include[] = "http://tefor.net/main/specific/nagging/nagging-menu.js";
    $include[] = "http://tefor.net/main/specific/style/specific.css";    
    
    $include[] = "http://tefor.net/main/style/newfont/personnalfont.css";
    $include[] = "http://tefor.net/main/specific/style/button.css";
    $include[] = "http://tefor.net/crispor/style/style.css";
    $include[] = "./style/style.css";
    
    foreach ($include as $inc_element)
    {
        if (pathinfo($inc_element,PATHINFO_EXTENSION)=="css")
            echo "<link rel='stylesheet' media='screen' type='text/css' href='$inc_element?v=2".filemtime($inc_element)."'/>
                 ";
        
        elseif (pathinfo($inc_element,PATHINFO_EXTENSION)=="js")
            echo "<script src='$inc_element?v=2".filemtime($inc_element)."'></script>
                 ";
    }
?>
