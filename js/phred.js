$(document).ready(function() {

    $('#input-phredScore').on('change keyup', function() {
        phredScore = Number($(this).val());

        errorProb = Math.pow(10, -phredScore / 10);
        baseAccuracy = 100 * (1 - errorProb)
        asciiValue = String.fromCharCode(phredScore + 33);
        
        $('#input-errorProbability').val(errorProb);
        $('#input-baseAccuracy').val(baseAccuracy);
        $('#input-asciiValue').val(asciiValue);
        
    });

    $('#input-errorProbability').on('change keyup', function() {
        errorProb = Number($(this).val());

        phredScore = -10 * Math.log10(errorProb);
        baseAccuracy = 100 * (1 - errorProb)
        asciiValue = String.fromCharCode(phredScore + 33);

        $('#input-phredScore').val(phredScore);
        $('#input-baseAccuracy').val(baseAccuracy);
        $('#input-asciiValue').val(asciiValue);
    });

    $('#input-baseAccuracy').on('change keyup', function() {
        baseAccuracy = Number($(this).val())

        errorProb = (100 - baseAccuracy) / 100;
        phredScore = -10 * Math.log10(errorProb);
        asciiValue = String.fromCharCode(phredScore + 33);

        $('#input-phredScore').val(phredScore);
        $('#input-errorProbability').val(errorProb);
        $('#input-asciiValue').val(asciiValue);
    });

    $('#input-asciiValue').on('change keyup', function() {
        phredScore = $(this).val().charCodeAt(0) - 33;
        errorProb = Math.pow(10, -phredScore / 10);
        baseAccuracy = 100 * (1 - errorProb)

        $('#input-phredScore').val(phredScore);
        $('#input-errorProbability').val(errorProb);
        $('#input-baseAccuracy').val(baseAccuracy);
    });

});
