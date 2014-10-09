<?php
    $arguments = "";
    foreach ($_GET as $key=>$value)
    {
        $arguments.=$key."=".$value."&";
    }
    
    $filename = 'CRISPOR-'.date("Y-m-d H:i:s").'.html'; // of course find the exact filename....        
    header('Pragma: public');
    header('Expires: 0');
    header('Cache-Control: must-revalidate, post-check=0, pre-check=0');
    header('Cache-Control: private', false); // required for certain browsers 
    header('Content-Type: application/html');

    header('Content-Disposition: attachment; filename="'. basename($filename) . '";');
    //header('Content-Transfer-Encoding: binary');
    //header('Content-Length: ' . filesize($filename));

    readfile("http://tefor.net/crispor/crispor.cgi?$arguments");


    exit;
?>