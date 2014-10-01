<?php
    $title = "CRISPOR - TEFOR&#39;s CRISPR ONLINE TOOL";
    $metadescription = "CRISPOR - Free guide-RNAs designer tool - Design your Cas9 guide-RNAs and control the off-targets in a growing list of genomes";
        
    echo "<title>$title</title>";
    echo "<meta name='description' content='$metadescription'/>";
    echo "<meta http-equiv='Content-Type' content='text/html; charset=utf-8' />";
    echo "<meta property='fb:admins' content='692090743' />";
    echo "<meta property='og:title' content='$title' />";
    echo "<meta property='og:description' content='$metadescription'/>";
    echo "<meta property='og:type' content='website' />";
    echo "<meta property='og:url' content='http://tefor.net/crispor/crispor.cgi' />";
    echo "<meta property='og:image' content='http://tefor.net/main/images/blog/blog/259/20140717172157_Sans_titre.png' />";
    
    echo "<script src='//code.jquery.com/jquery-1.11.0.min.js'></script>";   
    echo "<script src='//code.jquery.com/ui/1.11.1/jquery-ui.min.js'></script>";   
    echo "<script src='https://apis.google.com/js/platform.js' async defer></script>";

    
    $include[] = "http://tefor.net/main/specific/nagging/style.css";
//    $include[] = "http://tefor.net/main/specific/nagging/nagging-menu.js";
    $include[] = "http://tefor.net/main/specific/style/specific.css";    
    $include[] = "http://tefor.net/main/style/style_general.css";
    $include[] = "http://tefor.net/main/style/newfont/personnalfont.css";
    $include[] = "http://tefor.net/crispor/style/style.css";
    
    foreach ($include as $inc_element)
    {
        if (pathinfo($inc_element,PATHINFO_EXTENSION)=="css")
            echo "<link rel='stylesheet' media='screen' type='text/css' href='$inc_element?v=".filemtime($inc_element)."'/>";
        
        elseif (pathinfo($inc_element,PATHINFO_EXTENSION)=="js")
            echo "<script src='$inc_element?v=".filemtime($inc_element)."'></script>";        
    }
?>
