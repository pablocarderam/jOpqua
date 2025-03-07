using StaticArrays

# Model parameter setup initializers

function newPathogenType(
    id::String;
    template::PathogenType=DEFAULT_PATHOGEN_TYPE,
    num_loci::Union{Nothing,Int64}=nothing,
    possible_alleles::Union{Nothing,String}=nothing,

    # Each element takes seq argument, returns Float64
    mutantEstablishmentCoefficient::Union{Nothing,Function}=nothing,
    clearanceCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionCoefficient::Union{Nothing,Function}=nothing,
    recombinantEstablishmentCoefficient::Union{Nothing,Function}=nothing,
    contactCoefficient::Union{Nothing,Function}=nothing,
    responseLossCoefficient::Union{Nothing,Function}=nothing,
    birthCoefficient::Union{Nothing,Function}=nothing,
    deathCoefficient::Union{Nothing,Function}=nothing,
    transitionCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionCoefficient::Union{Nothing,Function}=nothing,
    receiveContactCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessCoefficient::Union{Nothing,Function}=nothing,
    mutationsUponInfectionCoefficient::Union{Nothing,Function}=nothing,
    recombinationsUponInfectionCoefficient::Union{Nothing,Function}=nothing,
    inoculumCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthCoefficient::Union{Nothing,Function}=nothing,
)

    isnothing(num_loci) ? num_loci = template.num_loci : num_loci = num_loci
    isnothing(possible_alleles) ? possible_alleles = template.possible_alleles : possible_alleles = possible_alleles

    isnothing(mutantEstablishmentCoefficient) ? mutantEstablishmentCoefficient = template.coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentCoefficient = mutantEstablishmentCoefficient
    isnothing(clearanceCoefficient) ? clearanceCoefficient = template.coefficient_functions[CLEARANCE] : clearanceCoefficient = clearanceCoefficient
    isnothing(responseAcquisitionCoefficient) ? responseAcquisitionCoefficient = template.coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionCoefficient = responseAcquisitionCoefficient
    isnothing(recombinantEstablishmentCoefficient) ? recombinantEstablishmentCoefficient = template.coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentCoefficient = recombinantEstablishmentCoefficient
    isnothing(contactCoefficient) ? contactCoefficient = template.coefficient_functions[CONTACT] : contactCoefficient = contactCoefficient
    isnothing(responseLossCoefficient) ? responseLossCoefficient = template.coefficient_functions[RESPONSE_LOSS] : responseLossCoefficient = responseLossCoefficient
    isnothing(birthCoefficient) ? birthCoefficient = template.coefficient_functions[BIRTH] : birthCoefficient = birthCoefficient
    isnothing(deathCoefficient) ? deathCoefficient = template.coefficient_functions[DEATH] : deathCoefficient = deathCoefficient
    isnothing(transitionCoefficient) ? transitionCoefficient = template.coefficient_functions[TRANSITION] : transitionCoefficient = transitionCoefficient
    isnothing(receiveTransitionCoefficient) ? receiveTransitionCoefficient = template.coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionCoefficient = receiveTransitionCoefficient
    isnothing(receiveContactCoefficient) ? receiveContactCoefficient = template.coefficient_functions[RECEIVE_CONTACT] : receiveContactCoefficient = receiveContactCoefficient
    isnothing(intrahostFitnessCoefficient) ? intrahostFitnessCoefficient = template.coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessCoefficient = intrahostFitnessCoefficient
    isnothing(mutationsUponInfectionCoefficient) ? mutationsUponInfectionCoefficient = template.coefficient_functions[MUTATIONS_UPON_INFECTION] : mutationsUponInfectionCoefficient = mutationsUponInfectionCoefficient
    isnothing(recombinationsUponInfectionCoefficient) ? recombinationsUponInfectionCoefficient = template.coefficient_functions[RECOMBINATIONS_UPON_INFECTION] : recombinationsUponInfectionCoefficient = recombinationsUponInfectionCoefficient
    isnothing(inoculumCoefficient) ? inoculumCoefficient = template.coefficient_functions[INOCULUM] : inoculumCoefficient = inoculumCoefficient
    isnothing(transmissionEfficiencyCoefficient) ? transmissionEfficiencyCoefficient = template.coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyCoefficient = transmissionEfficiencyCoefficient
    isnothing(verticalTransmissionCoefficient) ? verticalTransmissionCoefficient = template.coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionCoefficient = verticalTransmissionCoefficient
    isnothing(responseAcquisitionUponClearanceCoefficient) ? responseAcquisitionUponClearanceCoefficient = template.coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceCoefficient = responseAcquisitionUponClearanceCoefficient
    isnothing(responseInheritanceCoefficient) ? responseInheritanceCoefficient = template.coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceCoefficient = responseInheritanceCoefficient
    isnothing(hostMutationsUponBirthCoefficient) ? hostMutationsUponBirthCoefficient = template.coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthCoefficient = hostMutationsUponBirthCoefficient
    isnothing(hostRecombinationsUponBirthCoefficient) ? hostRecombinationsUponBirthCoefficient = template.coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthCoefficient = hostRecombinationsUponBirthCoefficient

    return PathogenType(
        id,
        num_loci,
        possible_alleles,
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentCoefficient, clearanceCoefficient, responseAcquisitionCoefficient,
            recombinantEstablishmentCoefficient, contactCoefficient, responseLossCoefficient,
            birthCoefficient, deathCoefficient, transitionCoefficient,
            receiveTransitionCoefficient, receiveContactCoefficient, intrahostFitnessCoefficient,
            mutationsUponInfectionCoefficient, recombinationsUponInfectionCoefficient, inoculumCoefficient,
            transmissionEfficiencyCoefficient, verticalTransmissionCoefficient, responseAcquisitionUponClearanceCoefficient,
            responseInheritanceCoefficient, hostMutationsUponBirthCoefficient, hostRecombinationsUponBirthCoefficient,
        ],
    )
end

function newResponseType(
    id::String;
    template::ResponseType=DEFAULT_RESPONSE_TYPE,
    reactivityCoefficient::Union{Nothing,Function}=nothing,
    # takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient

    # Each takes host, imprinted, matured sequences and returns Float64 coefficient
    mutantEstablishmentStaticCoefficient::Union{Nothing,Function}=nothing,
    clearanceStaticCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionStaticCoefficient::Union{Nothing,Function}=nothing,
    recombinantEstablishmentStaticCoefficient::Union{Nothing,Function}=nothing,
    contactStaticCoefficient::Union{Nothing,Function}=nothing,
    responseLossStaticCoefficient::Union{Nothing,Function}=nothing,
    birthStaticCoefficient::Union{Nothing,Function}=nothing,
    deathStaticCoefficient::Union{Nothing,Function}=nothing,
    transitionStaticCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionStaticCoefficient::Union{Nothing,Function}=nothing,
    receiveContactStaticCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessStaticCoefficient::Union{Nothing,Function}=nothing,
    mutationsUponInfectionStaticCoefficient::Union{Nothing,Function}=nothing,
    recombinationsUponInfectionStaticCoefficient::Union{Nothing,Function}=nothing,
    inoculumStaticCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyStaticCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionStaticCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceStaticCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceStaticCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthStaticCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthStaticCoefficient::Union{Nothing,Function}=nothing,

    # Each takes host, imprinted, matured, and infecting sequences and returns Float64 coefficient
    mutantEstablishmentSpecificCoefficient::Union{Nothing,Function}=nothing,
    clearanceSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionSpecificCoefficient::Union{Nothing,Function}=nothing,
    recombinantEstablishmentSpecificCoefficient::Union{Nothing,Function}=nothing,
    contactSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseLossSpecificCoefficient::Union{Nothing,Function}=nothing,
    birthSpecificCoefficient::Union{Nothing,Function}=nothing,
    deathSpecificCoefficient::Union{Nothing,Function}=nothing,
    transitionSpecificCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionSpecificCoefficient::Union{Nothing,Function}=nothing,
    receiveContactSpecificCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessSpecificCoefficient::Union{Nothing,Function}=nothing,
    mutationsUponInfectionSpecificCoefficient::Union{Nothing,Function}=nothing,
    recombinationsUponInfectionSpecificCoefficient::Union{Nothing,Function}=nothing,
    inoculumSpecificCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencySpecificCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceSpecificCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceSpecificCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthSpecificCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthSpecificCoefficient::Union{Nothing,Function}=nothing,
)

    isnothing(reactivityCoefficient) ? reactivityCoefficient = template.reactivityCoefficient : reactivityCoefficient = reactivityCoefficient

    isnothing(mutantEstablishmentStaticCoefficient) ? mutantEstablishmentStaticCoefficient = template.static_coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentStaticCoefficient = mutantEstablishmentStaticCoefficient
    isnothing(clearanceStaticCoefficient) ? clearanceStaticCoefficient = template.static_coefficient_functions[CLEARANCE] : clearanceStaticCoefficient = clearanceStaticCoefficient
    isnothing(responseAcquisitionStaticCoefficient) ? responseAcquisitionStaticCoefficient = template.static_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionStaticCoefficient = responseAcquisitionStaticCoefficient
    isnothing(recombinantEstablishmentStaticCoefficient) ? recombinantEstablishmentStaticCoefficient = template.static_coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentStaticCoefficient = recombinantEstablishmentStaticCoefficient
    isnothing(contactStaticCoefficient) ? contactStaticCoefficient = template.static_coefficient_functions[CONTACT] : contactStaticCoefficient = contactStaticCoefficient
    isnothing(responseLossStaticCoefficient) ? responseLossStaticCoefficient = template.static_coefficient_functions[RESPONSE_LOSS] : responseLossStaticCoefficient = responseLossStaticCoefficient
    isnothing(birthStaticCoefficient) ? birthStaticCoefficient = template.static_coefficient_functions[BIRTH] : birthStaticCoefficient = birthStaticCoefficient
    isnothing(deathStaticCoefficient) ? deathStaticCoefficient = template.static_coefficient_functions[DEATH] : deathStaticCoefficient = deathStaticCoefficient
    isnothing(transitionStaticCoefficient) ? transitionStaticCoefficient = template.static_coefficient_functions[TRANSITION] : transitionStaticCoefficient = transitionStaticCoefficient
    isnothing(receiveTransitionStaticCoefficient) ? receiveTransitionStaticCoefficient = template.static_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionStaticCoefficient = receiveTransitionStaticCoefficient
    isnothing(receiveContactStaticCoefficient) ? receiveContactStaticCoefficient = template.static_coefficient_functions[RECEIVE_CONTACT] : receiveContactStaticCoefficient = receiveContactStaticCoefficient
    isnothing(intrahostFitnessStaticCoefficient) ? intrahostFitnessStaticCoefficient = template.static_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessStaticCoefficient = intrahostFitnessStaticCoefficient
    isnothing(mutationsUponInfectionStaticCoefficient) ? mutationsUponInfectionStaticCoefficient = template.static_coefficient_functions[MUTATIONS_UPON_INFECTION] : mutationsUponInfectionStaticCoefficient = mutationsUponInfectionStaticCoefficient
    isnothing(recombinationsUponInfectionStaticCoefficient) ? recombinationsUponInfectionStaticCoefficient = template.static_coefficient_functions[RECOMBINATIONS_UPON_INFECTION] : recombinationsUponInfectionStaticCoefficient = recombinationsUponInfectionStaticCoefficient
    isnothing(inoculumStaticCoefficient) ? inoculumStaticCoefficient = template.static_coefficient_functions[INOCULUM] : inoculumStaticCoefficient = inoculumStaticCoefficient
    isnothing(transmissionEfficiencyStaticCoefficient) ? transmissionEfficiencyStaticCoefficient = template.static_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyStaticCoefficient = transmissionEfficiencyStaticCoefficient
    isnothing(verticalTransmissionStaticCoefficient) ? verticalTransmissionStaticCoefficient = template.static_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionStaticCoefficient = verticalTransmissionStaticCoefficient
    isnothing(responseAcquisitionUponClearanceStaticCoefficient) ? responseAcquisitionUponClearanceStaticCoefficient = template.static_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceStaticCoefficient = responseAcquisitionUponClearanceStaticCoefficient
    isnothing(responseInheritanceStaticCoefficient) ? responseInheritanceStaticCoefficient = template.static_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceStaticCoefficient = responseInheritanceStaticCoefficient
    isnothing(hostMutationsUponBirthStaticCoefficient) ? hostMutationsUponBirthStaticCoefficient = template.static_coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthStaticCoefficient = hostMutationsUponBirthStaticCoefficient
    isnothing(hostRecombinationsUponBirthStaticCoefficient) ? hostRecombinationsUponBirthStaticCoefficient = template.static_coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthStaticCoefficient = hostRecombinationsUponBirthStaticCoefficient

    isnothing(mutantEstablishmentSpecificCoefficient) ? mutantEstablishmentSpecificCoefficient = template.specific_coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentSpecificCoefficient = mutantEstablishmentSpecificCoefficient
    isnothing(clearanceSpecificCoefficient) ? clearanceSpecificCoefficient = template.specific_coefficient_functions[CLEARANCE] : clearanceSpecificCoefficient = clearanceSpecificCoefficient
    isnothing(responseAcquisitionSpecificCoefficient) ? responseAcquisitionSpecificCoefficient = template.specific_coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionSpecificCoefficient = responseAcquisitionSpecificCoefficient
    isnothing(recombinantEstablishmentSpecificCoefficient) ? recombinantEstablishmentSpecificCoefficient = template.specific_coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentSpecificCoefficient = recombinantEstablishmentSpecificCoefficient
    isnothing(contactSpecificCoefficient) ? contactSpecificCoefficient = template.specific_coefficient_functions[CONTACT] : contactSpecificCoefficient = contactSpecificCoefficient
    isnothing(responseLossSpecificCoefficient) ? responseLossSpecificCoefficient = template.specific_coefficient_functions[RESPONSE_LOSS] : responseLossSpecificCoefficient = responseLossSpecificCoefficient
    isnothing(birthSpecificCoefficient) ? birthSpecificCoefficient = template.specific_coefficient_functions[BIRTH] : birthSpecificCoefficient = birthSpecificCoefficient
    isnothing(deathSpecificCoefficient) ? deathSpecificCoefficient = template.specific_coefficient_functions[DEATH] : deathSpecificCoefficient = deathSpecificCoefficient
    isnothing(transitionSpecificCoefficient) ? transitionSpecificCoefficient = template.specific_coefficient_functions[TRANSITION] : transitionSpecificCoefficient = transitionSpecificCoefficient
    isnothing(receiveTransitionSpecificCoefficient) ? receiveTransitionSpecificCoefficient = template.specific_coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionSpecificCoefficient = receiveTransitionSpecificCoefficient
    isnothing(receiveContactSpecificCoefficient) ? receiveContactSpecificCoefficient = template.specific_coefficient_functions[RECEIVE_CONTACT] : receiveContactSpecificCoefficient = receiveContactSpecificCoefficient
    isnothing(intrahostFitnessSpecificCoefficient) ? intrahostFitnessSpecificCoefficient = template.specific_coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessSpecificCoefficient = intrahostFitnessSpecificCoefficient
    isnothing(mutationsUponInfectionSpecificCoefficient) ? mutationsUponInfectionSpecificCoefficient = template.specific_coefficient_functions[MUTATIONS_UPON_INFECTION] : mutationsUponInfectionSpecificCoefficient = mutationsUponInfectionSpecificCoefficient
    isnothing(recombinationsUponInfectionSpecificCoefficient) ? recombinationsUponInfectionSpecificCoefficient = template.specific_coefficient_functions[RECOMBINATIONS_UPON_INFECTION] : recombinationsUponInfectionSpecificCoefficient = recombinationsUponInfectionSpecificCoefficient
    isnothing(inoculumSpecificCoefficient) ? inoculumSpecificCoefficient = template.specific_coefficient_functions[INOCULUM] : inoculumSpecificCoefficient = inoculumSpecificCoefficient
    isnothing(transmissionEfficiencySpecificCoefficient) ? transmissionEfficiencySpecificCoefficient = template.specific_coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencySpecificCoefficient = transmissionEfficiencySpecificCoefficient
    isnothing(verticalTransmissionSpecificCoefficient) ? verticalTransmissionSpecificCoefficient = template.specific_coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionSpecificCoefficient = verticalTransmissionSpecificCoefficient
    isnothing(responseAcquisitionUponClearanceSpecificCoefficient) ? responseAcquisitionUponClearanceSpecificCoefficient = template.specific_coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceSpecificCoefficient = responseAcquisitionUponClearanceSpecificCoefficient
    isnothing(responseInheritanceSpecificCoefficient) ? responseInheritanceSpecificCoefficient = template.specific_coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceSpecificCoefficient = responseInheritanceSpecificCoefficient
    isnothing(hostMutationsUponBirthSpecificCoefficient) ? hostMutationsUponBirthSpecificCoefficient = template.specific_coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthSpecificCoefficient = hostMutationsUponBirthSpecificCoefficient
    isnothing(hostRecombinationsUponBirthSpecificCoefficient) ? hostRecombinationsUponBirthSpecificCoefficient = template.specific_coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthSpecificCoefficient = hostRecombinationsUponBirthSpecificCoefficient

    return ResponseType(
        id,
        reactivityCoefficient,
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentStaticCoefficient, clearanceStaticCoefficient, responseAcquisitionStaticCoefficient,
            recombinantEstablishmentStaticCoefficient, contactStaticCoefficient, responseLossStaticCoefficient,
            birthStaticCoefficient, deathStaticCoefficient, transitionStaticCoefficient,
            receiveTransitionStaticCoefficient, receiveContactStaticCoefficient, intrahostFitnessStaticCoefficient,
            mutationsUponInfectionStaticCoefficient, recombinationsUponInfectionStaticCoefficient, inoculumStaticCoefficient,
            transmissionEfficiencyStaticCoefficient, verticalTransmissionStaticCoefficient, responseAcquisitionUponClearanceStaticCoefficient,
            responseInheritanceStaticCoefficient, hostMutationsUponBirthStaticCoefficient, hostRecombinationsUponBirthStaticCoefficient,
        ],
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentSpecificCoefficient, clearanceSpecificCoefficient, responseAcquisitionSpecificCoefficient,
            recombinantEstablishmentSpecificCoefficient, contactSpecificCoefficient, responseLossSpecificCoefficient,
            birthSpecificCoefficient, deathSpecificCoefficient, transitionSpecificCoefficient,
            receiveTransitionSpecificCoefficient, receiveContactSpecificCoefficient, intrahostFitnessSpecificCoefficient,
            mutationsUponInfectionSpecificCoefficient, recombinationsUponInfectionSpecificCoefficient, inoculumSpecificCoefficient,
            transmissionEfficiencySpecificCoefficient, verticalTransmissionSpecificCoefficient, responseAcquisitionUponClearanceSpecificCoefficient,
            responseInheritanceSpecificCoefficient, hostMutationsUponBirthSpecificCoefficient, hostRecombinationsUponBirthSpecificCoefficient,
        ],
    )
end

function newHostType(
    id::String;
    template::HostType=DEFAULT_HOST_TYPE,
    num_loci::Union{Nothing,Int64}=nothing,
    possible_alleles::Union{Nothing,String}=nothing,

    # Each element takes seq argument, returns Float64
    mutantEstablishmentCoefficient::Union{Nothing,Function}=nothing,
    clearanceCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionCoefficient::Union{Nothing,Function}=nothing,
    recombinantEstablishmentCoefficient::Union{Nothing,Function}=nothing,
    contactCoefficient::Union{Nothing,Function}=nothing,
    responseLossCoefficient::Union{Nothing,Function}=nothing,
    birthCoefficient::Union{Nothing,Function}=nothing,
    deathCoefficient::Union{Nothing,Function}=nothing,
    transitionCoefficient::Union{Nothing,Function}=nothing,
    receiveTransitionCoefficient::Union{Nothing,Function}=nothing,
    receiveContactCoefficient::Union{Nothing,Function}=nothing,
    intrahostFitnessCoefficient::Union{Nothing,Function}=nothing,
    mutationsUponInfectionCoefficient::Union{Nothing,Function}=nothing,
    recombinationsUponInfectionCoefficient::Union{Nothing,Function}=nothing,
    inoculumCoefficient::Union{Nothing,Function}=nothing,
    transmissionEfficiencyCoefficient::Union{Nothing,Function}=nothing,
    verticalTransmissionCoefficient::Union{Nothing,Function}=nothing,
    responseAcquisitionUponClearanceCoefficient::Union{Nothing,Function}=nothing,
    responseInheritanceCoefficient::Union{Nothing,Function}=nothing,
    hostMutationsUponBirthCoefficient::Union{Nothing,Function}=nothing,
    hostRecombinationsUponBirthCoefficient::Union{Nothing,Function}=nothing,
)

    isnothing(num_loci) ? num_loci = template.num_loci : num_loci = num_loci
    isnothing(possible_alleles) ? possible_alleles = template.possible_alleles : possible_alleles = possible_alleles

    isnothing(mutantEstablishmentCoefficient) ? mutantEstablishmentCoefficient = template.coefficient_functions[MUTANT_ESTABLISHMENT] : mutantEstablishmentCoefficient = mutantEstablishmentCoefficient
    isnothing(clearanceCoefficient) ? clearanceCoefficient = template.coefficient_functions[CLEARANCE] : clearanceCoefficient = clearanceCoefficient
    isnothing(responseAcquisitionCoefficient) ? responseAcquisitionCoefficient = template.coefficient_functions[RESPONSE_ACQUISITION] : responseAcquisitionCoefficient = responseAcquisitionCoefficient
    isnothing(recombinantEstablishmentCoefficient) ? recombinantEstablishmentCoefficient = template.coefficient_functions[RECOMBINANT_ESTABLISHMENT] : recombinantEstablishmentCoefficient = recombinantEstablishmentCoefficient
    isnothing(contactCoefficient) ? contactCoefficient = template.coefficient_functions[CONTACT] : contactCoefficient = contactCoefficient
    isnothing(responseLossCoefficient) ? responseLossCoefficient = template.coefficient_functions[RESPONSE_LOSS] : responseLossCoefficient = responseLossCoefficient
    isnothing(birthCoefficient) ? birthCoefficient = template.coefficient_functions[BIRTH] : birthCoefficient = birthCoefficient
    isnothing(deathCoefficient) ? deathCoefficient = template.coefficient_functions[DEATH] : deathCoefficient = deathCoefficient
    isnothing(transitionCoefficient) ? transitionCoefficient = template.coefficient_functions[TRANSITION] : transitionCoefficient = transitionCoefficient
    isnothing(receiveTransitionCoefficient) ? receiveTransitionCoefficient = template.coefficient_functions[RECEIVE_TRANSITION] : receiveTransitionCoefficient = receiveTransitionCoefficient
    isnothing(receiveContactCoefficient) ? receiveContactCoefficient = template.coefficient_functions[RECEIVE_CONTACT] : receiveContactCoefficient = receiveContactCoefficient
    isnothing(intrahostFitnessCoefficient) ? intrahostFitnessCoefficient = template.coefficient_functions[INTRAHOST_FITNESS] : intrahostFitnessCoefficient = intrahostFitnessCoefficient
    isnothing(mutationsUponInfectionCoefficient) ? mutationsUponInfectionCoefficient = template.coefficient_functions[MUTATIONS_UPON_INFECTION] : mutationsUponInfectionCoefficient = mutationsUponInfectionCoefficient
    isnothing(recombinationsUponInfectionCoefficient) ? recombinationsUponInfectionCoefficient = template.coefficient_functions[RECOMBINATIONS_UPON_INFECTION] : recombinationsUponInfectionCoefficient = recombinationsUponInfectionCoefficient
    isnothing(inoculumCoefficient) ? inoculumCoefficient = template.coefficient_functions[INOCULUM] : inoculumCoefficient = inoculumCoefficient
    isnothing(transmissionEfficiencyCoefficient) ? transmissionEfficiencyCoefficient = template.coefficient_functions[TRANSMISSION_EFFICIENCY] : transmissionEfficiencyCoefficient = transmissionEfficiencyCoefficient
    isnothing(verticalTransmissionCoefficient) ? verticalTransmissionCoefficient = template.coefficient_functions[VERTICAL_TRANSMISSION] : verticalTransmissionCoefficient = verticalTransmissionCoefficient
    isnothing(responseAcquisitionUponClearanceCoefficient) ? responseAcquisitionUponClearanceCoefficient = template.coefficient_functions[RESPONSE_ACQUISITION_UPON_CLEARANCE] : responseAcquisitionUponClearanceCoefficient = responseAcquisitionUponClearanceCoefficient
    isnothing(responseInheritanceCoefficient) ? responseInheritanceCoefficient = template.coefficient_functions[RESPONSE_INHERITANCE] : responseInheritanceCoefficient = responseInheritanceCoefficient
    isnothing(hostMutationsUponBirthCoefficient) ? hostMutationsUponBirthCoefficient = template.coefficient_functions[HOST_MUTATIONS_UPON_BIRTH] : hostMutationsUponBirthCoefficient = hostMutationsUponBirthCoefficient
    isnothing(hostRecombinationsUponBirthCoefficient) ? hostRecombinationsUponBirthCoefficient = template.coefficient_functions[HOST_RECOMBINATIONS_UPON_BIRTH] : hostRecombinationsUponBirthCoefficient = hostRecombinationsUponBirthCoefficient

    return HostType(
        id,
        num_loci,
        possible_alleles,
        SA[ # order defined in COEFFICIENTS
            mutantEstablishmentCoefficient, clearanceCoefficient, responseAcquisitionCoefficient,
            recombinantEstablishmentCoefficient, contactCoefficient, responseLossCoefficient,
            birthCoefficient, deathCoefficient, transitionCoefficient,
            receiveTransitionCoefficient, receiveContactCoefficient, intrahostFitnessCoefficient,
            mutationsUponInfectionCoefficient, recombinationsUponInfectionCoefficient, inoculumCoefficient,
            transmissionEfficiencyCoefficient, verticalTransmissionCoefficient, responseAcquisitionUponClearanceCoefficient,
            responseInheritanceCoefficient, hostMutationsUponBirthCoefficient, hostRecombinationsUponBirthCoefficient,
        ],
    )
end

function newPopulationType(
    id::String;
    template::PopulationType=DEFAULT_POPULATION_TYPE,
    constant_contact_density::Union{Nothing,Bool}=nothing,
    constant_transition_density::Union{Nothing,Bool}=nothing,
    host_sexual_reproduction::Union{Nothing,Bool}=nothing,
    hostSexualCompatibility::Union{Nothing,Function}=nothing,

    # Rate coefficients:
    mutant_establishment_coefficient::Union{Nothing,Float64}=nothing,
    clearance_coefficient::Union{Nothing,Float64}=nothing,
    response_acquisition_coefficient::Union{Nothing,Float64}=nothing,
    recombinant_establishment_coefficient::Union{Nothing,Float64}=nothing,
    contact_coefficient::Union{Nothing,Float64}=nothing,
    response_loss_coefficient::Union{Nothing,Float64}=nothing,
    birth_coefficient::Union{Nothing,Float64}=nothing,
    death_coefficient::Union{Nothing,Float64}=nothing,
    transition_coefficient::Union{Nothing,Float64}=nothing,
    intrahost_fitness_coefficient::Union{Nothing,Float64}=nothing,
    mutations_upon_infection_coefficient::Union{Nothing,Float64}=nothing,
    recombinations_upon_infection_coefficient::Union{Nothing,Float64}=nothing,
    inoculum_coefficient::Union{Nothing,Float64}=nothing,
    transmission_efficiency_coefficient::Union{Nothing,Float64}=nothing,
    vertical_transmission_coefficient::Union{Nothing,Float64}=nothing,
    response_acquisition_upon_clearance_coefficient::Union{Nothing,Float64}=nothing,
    response_inheritance_coefficient::Union{Nothing,Float64}=nothing,
    host_mutations_upon_birth_coefficient::Union{Nothing,Float64}=nothing,
    host_recombinations_upon_birth_coefficient::Union{Nothing,Float64}=nothing,

    # Receive coefficients
    receive_transition_coefficient::Union{Nothing,Float64}=nothing,
    receive_contact_coefficient::Union{Nothing,Float64}=nothing,

    # Other functions
    pathogenFractions::Union{Nothing,Function}=nothing,
    # Takes Host entity and Population's weightedInteraction function,
    # returns vector with fractional representation of each pathogen present
    weightedInteraction::Union{Nothing,Function}=nothing,
    # Takes Pathogen entity, Host entity, and event number;
    # returns aggregated response coefficient against that Pathogen for that event
    transmissionEfficiency::Union{Nothing,Function}=nothing,
    # Takes Pathogen and Host entities,
    # returns probability that a contact results in successful infection given the Responses in Host

    developResponses::Union{Nothing,Function}=nothing,
    # takes in Pathogen, Host, population's Dict of Responses, population type's dictionary of ResponseTypes, and birth time as arguments, returns Response entities to be added
    # (this handles how many and which responses to choose when adding a response to a host)

    response_types::Union{Nothing,Dict{String,ResponseType}}=nothing,
)

    isnothing(constant_contact_density) ? constant_contact_density = template.constant_contact_density : constant_contact_density = constant_contact_density
    isnothing(constant_transition_density) ? constant_transition_density = template.constant_transition_density : constant_transition_density = constant_transition_density
    isnothing(host_sexual_reproduction) ? host_sexual_reproduction = template.host_sexual_reproduction : host_sexual_reproduction = host_sexual_reproduction

    isnothing(hostSexualCompatibility) ? hostSexualCompatibility = template.hostSexualCompatibility : hostSexualCompatibility = hostSexualCompatibility

    isnothing(mutant_establishment_coefficient) ? mutant_establishment_coefficient = template.base_coefficients[MUTANT_ESTABLISHMENT] : mutant_establishment_coefficient = mutant_establishment_coefficient
    isnothing(clearance_coefficient) ? clearance_coefficient = template.base_coefficients[CLEARANCE] : clearance_coefficient = clearance_coefficient
    isnothing(response_acquisition_coefficient) ? response_acquisition_coefficient = template.base_coefficients[RESPONSE_ACQUISITION] : response_acquisition_coefficient = response_acquisition_coefficient
    isnothing(recombinant_establishment_coefficient) ? recombinant_establishment_coefficient = template.base_coefficients[RECOMBINANT_ESTABLISHMENT] : recombinant_establishment_coefficient = recombinant_establishment_coefficient
    isnothing(contact_coefficient) ? contact_coefficient = template.base_coefficients[CONTACT] : contact_coefficient = contact_coefficient
    isnothing(response_loss_coefficient) ? response_loss_coefficient = template.base_coefficients[RESPONSE_LOSS] : response_loss_coefficient = response_loss_coefficient
    isnothing(birth_coefficient) ? birth_coefficient = template.base_coefficients[BIRTH] : birth_coefficient = birth_coefficient
    isnothing(death_coefficient) ? death_coefficient = template.base_coefficients[DEATH] : death_coefficient = death_coefficient
    isnothing(transition_coefficient) ? transition_coefficient = template.base_coefficients[TRANSITION] : transition_coefficient = transition_coefficient
    isnothing(receive_transition_coefficient) ? receive_transition_coefficient = template.base_coefficients[RECEIVE_TRANSITION] : receive_transition_coefficient = receive_transition_coefficient
    isnothing(receive_contact_coefficient) ? receive_contact_coefficient = template.base_coefficients[RECEIVE_CONTACT] : receive_contact_coefficient = receive_contact_coefficient
    isnothing(intrahost_fitness_coefficient) ? intrahost_fitness_coefficient = template.base_coefficients[INTRAHOST_FITNESS] : intrahost_fitness_coefficient = intrahost_fitness_coefficient
    isnothing(mutations_upon_infection_coefficient) ? mutations_upon_infection_coefficient = template.base_coefficients[MUTATIONS_UPON_INFECTION] : mutations_upon_infection_coefficient = mutations_upon_infection_coefficient
    isnothing(recombinations_upon_infection_coefficient) ? recombinations_upon_infection_coefficient = template.base_coefficients[RECOMBINATIONS_UPON_INFECTION] : recombinations_upon_infection_coefficient = recombinations_upon_infection_coefficient
    isnothing(inoculum_coefficient) ? inoculum_coefficient = template.base_coefficients[INOCULUM] : inoculum_coefficient = inoculum_coefficient
    isnothing(transmission_efficiency_coefficient) ? transmission_efficiency_coefficient = template.base_coefficients[TRANSMISSION_EFFICIENCY] : transmission_efficiency_coefficient = transmission_efficiency_coefficient
    isnothing(vertical_transmission_coefficient) ? vertical_transmission_coefficient = template.base_coefficients[VERTICAL_TRANSMISSION] : vertical_transmission_coefficient = vertical_transmission_coefficient
    isnothing(response_acquisition_upon_clearance_coefficient) ? response_acquisition_upon_clearance_coefficient = template.base_coefficients[RESPONSE_ACQUISITION_UPON_CLEARANCE] : response_acquisition_upon_clearance_coefficient = response_acquisition_upon_clearance_coefficient
    isnothing(response_inheritance_coefficient) ? response_inheritance_coefficient = template.base_coefficients[RESPONSE_INHERITANCE] : response_inheritance_coefficient = response_inheritance_coefficient
    isnothing(host_mutations_upon_birth_coefficient) ? host_mutations_upon_birth_coefficient = template.base_coefficients[HOST_MUTATIONS_UPON_BIRTH] : host_mutations_upon_birth_coefficient = host_mutations_upon_birth_coefficient
    isnothing(host_recombinations_upon_birth_coefficient) ? host_recombinations_upon_birth_coefficient = template.base_coefficients[HOST_RECOMBINATIONS_UPON_BIRTH] : host_recombinations_upon_birth_coefficient = host_recombinations_upon_birth_coefficient

    isnothing(pathogenFractions) ? pathogenFractions = template.pathogenFractions : pathogenFractions = pathogenFractions
    isnothing(weightedInteraction) ? weightedInteraction = template.weightedInteraction : weightedInteraction = weightedInteraction
    isnothing(transmissionEfficiency) ? transmissionEfficiency = template.transmissionEfficiency : transmissionEfficiency = transmissionEfficiency
    isnothing(developResponses) ? developResponses = template.developResponses : developResponses = developResponses

    isnothing(response_types) ? response_types = template.response_types : response_types = response_types

    return PopulationType(
        id,
        constant_contact_density,
        constant_transition_density,
        host_sexual_reproduction,
        hostSexualCompatibility,
        SA[ # order defined in COEFFICIENTS
            mutant_establishment_coefficient, clearance_coefficient, response_acquisition_coefficient,
            recombinant_establishment_coefficient, contact_coefficient, response_loss_coefficient,
            birth_coefficient, death_coefficient, transition_coefficient,
            receive_transition_coefficient, receive_contact_coefficient, intrahost_fitness_coefficient,
            mutations_upon_infection_coefficient, recombinations_upon_infection_coefficient, inoculum_coefficient,
            transmission_efficiency_coefficient, vertical_transmission_coefficient, response_acquisition_upon_clearance_coefficient,
            response_inheritance_coefficient, host_mutations_upon_birth_coefficient, host_recombinations_upon_birth_coefficient,
        ],
        pathogenFractions,
        weightedInteraction,
        transmissionEfficiency,
        developResponses,
        response_types,
    )
end
