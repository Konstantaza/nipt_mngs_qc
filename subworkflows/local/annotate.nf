// =============================================================================
// ANNOTATE — subworkflow таксономической аннотации (интерфейс v1.4.0)
//
//   ANNOTATE(reads, dbsheet) -> profiles: (meta, profile.tsv,
//                                          profile.json, native)
//
//   reads    — [meta(id, fraction, single_end), r1, r2, se]
//   dbsheet  — канал строк [db_id, annotator, path, role] из --dbsheet
//
// meta на выходе дополняется db_id, role, annotator. Каждый выход
// адаптера проходит VALIDATE_PROFILE (контракт profile v2,
// docs/profile_contract.md).
//
// РЕЕСТР аннотаторов — ветки branch ниже. Добавление аннотатора:
//   1) modules/local/annotate_<id>.nf (процесс с emit: annotated той же
//      формы: (meta, profile.tsv, profile.json, native));
//   2) bin/<id>_to_profile.py (конвертер в контракт profile v2);
//   3) строка в dbsheet (annotator=<id>);
//   4) одна ветка в branch ниже + include процесса.
// main.nf и потребители (COMPARE, QC_GATE) при этом НЕ меняются.
// =============================================================================

include { ANNOTATE_KRAKEN2_BRACKEN } from '../../modules/local/annotate_kraken2_bracken'
include { ANNOTATE_STUB }            from '../../modules/local/annotate_stub'
include { ANNOTATE_METAX }           from '../../modules/local/annotate_metax'
include { VALIDATE_PROFILE }         from '../../modules/local/validate_profile'

workflow ANNOTATE {
    take:
    reads    // [meta, r1, r2, se]
    dbsheet  // канал maps [db_id, annotator, path, role]

    main:
    // fork: каждая фракция x каждая строка dbsheet; meta дополняется
    // полями БД. path уже резолвлен в file() при парсинге dbsheet.
    ch_in = reads
        .combine(dbsheet)
        .map { meta, r1, r2, se, dbrow ->
            [meta + [db_id: dbrow.db_id, role: dbrow.role,
                     annotator: dbrow.annotator],
             r1, r2, se, dbrow.path]
        }

    // Реестр аннотаторов: ветка по колонке annotator
    ch_branched = ch_in.branch { meta, r1, r2, se, dbpath ->
        kraken2_bracken: meta.annotator == 'kraken2_bracken'
        stub:            meta.annotator == 'stub'
        metax:           meta.annotator == 'metax'
        other:           true
    }
    // Неизвестный annotator — явная ошибка, а не молчаливая потеря
    // профилей образца (иначе сработало бы только FAIL/PIPELINE_ERROR
    // в QC gate без указания причины)
    ch_branched.other.subscribe { meta, r1, r2, se, dbpath ->
        throw new IllegalArgumentException(
            "ANNOTATE: неизвестный annotator '${meta.annotator}' " +
            "(db_id='${meta.db_id}'). Реестр: kraken2_bracken, stub, metax")
    }

    ANNOTATE_KRAKEN2_BRACKEN(ch_branched.kraken2_bracken)
    ANNOTATE_STUB(ch_branched.stub)
    ANNOTATE_METAX(ch_branched.metax)

    ch_annotated = ANNOTATE_KRAKEN2_BRACKEN.out.annotated
        .mix(ANNOTATE_STUB.out.annotated)
        .mix(ANNOTATE_METAX.out.annotated)

    VALIDATE_PROFILE(ch_annotated)

    emit:
    // (meta[id, fraction, single_end, db_id, role, annotator],
    //  profile.tsv, profile.json, native)
    profiles = VALIDATE_PROFILE.out.profiles
}
